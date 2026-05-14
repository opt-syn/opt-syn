classdef (Abstract) opt_manager_interface < handle
    %OPT_MANAGER_INTERFACE interface for the analysis and synthesis of
    %optimization/fixed point algorithms

    %
    %inheritance:
    % opt_analysis  < opt_manager_interface 
    % opt_synthesis < opt_manager_interface 
    
    properties
        sys;  %system (opt_system type)
        cons = [];
        vars = {}; %specifications
        iqc_op = {};
        specs = {};        

        %task: 'analysis' or 'synthesis': determined by subclass
        
        
        config = [];

        lmi = [];
        %other options
        task = 'generic';
    end
    
    methods
        function obj = opt_manager_interface(sys, config)
            %OPT_MANAGER_INTERFACE

            if nargin < 2
                config = opt_config();
            end

            obj.sys = sys;
            nop = length(obj.sys.op);
            obj.iqc_op = cell(nop, 1);
            obj.vars = struct('op', []);
            obj.vars.op =cell(nop, 1);   

            obj.config = config;
            

        end

        function lmi_handler = select_lmi(obj, sys)
            %SELECT_LMI select the lmi routines based on the system type
            %
            %Input:
            %   sys:    type of system (e.g. opt_system_switched)
            %   routine:    'analysis' or 'synthesis'
            
            tp = sys.get_type();
            clname = ['lmi_', obj.task, '_', tp];
            lmi_hand = str2func(clname);

            lmi_handler = lmi_hand(sys, obj.config);

        end



        function [vars, cons, objective, alg_psi, rho] = build_program(obj, specs)
            %BUILD_PROGRAM set up the algorithm analysis or synthesis problem

            if nargin < 2
                specs = obj.specs;
            end

            cons = obj.cons;
            vars = obj.vars;

            % [diss, cons] = obj.build_plant(cons);
            %alg_loop: used for debugging. The algorithm after signal 
            % transformationsbefore, but before cascade by the filters    

            [rho, sperf] = obj.perf_specs(specs);

            [iqc_data] = obj.iqc_op_all();
            [alg_psi, iqc_op, alg_loop] = obj.sys.build_plant(iqc_data, rho);



            [vars, cons] = obj.lmi.create_vars(vars, cons, alg_psi, sperf);

            %the dissipation can change
            [vars, cons, objective] = obj.cons_dynamic(vars, cons, alg_psi, iqc_op, specs);

        end

        function verdict = LMILAB(obj)
            %is LMILAB used?
            verdict = obj.config.LMILAB();
        end

        %% acquire solutions
        function [sol] = run(obj,  vars, cons, objective)
            %RUN: run the program
            sol = struct;
            
            
            if obj.LMILAB
                if ~isnumeric(objective)
                    lmis(cons, objective, 'c');
                end

                [lmi_out,info_out]=lmisolve(cons);
                
                
                sol.dia = lmi_out.dia;
                STATUS = (lmi_out.status || (sol.dia > obj.config.tol.dia));
                sol.info = info_out;


                %explicitly forming the blocks in recovery takes a long
                %time. Why is this the case?

                
                % ncons = length(cons.lmim);
                % sol.blocks = cell(ncons, 1);
                % sol.eb = zeros(ncons, 1);
                % for i = 1:length(cons.lmim)
                %     sol.blocks{i} = -double(double(cons.lmim(i), lmi_out));
                %     sol.eb(i) = min(eig(sol.blocks{i}));
                % end
                
            else %YALMIP
                opt = sdpsettings('verbose', obj.gen.config.verbose, 'solver', obj.config.gen.solver);
            
                t = optimize(lmi, objective, opt);
                
                STATUS = t.problem;
            end
            sol.status = STATUS;
            sol.lmi_out = lmi_out;

            
        end


        function obj = purge(obj);
            %PURGE: get rid of problem description
            obj.cons = [];
            obj.vars = {};
            obj.iqc = {};
            obj.specs = {};
        end
        
        function obj = add_specifications(obj, varargin)
            %concatenate the new performance specifications
            if iscell(varargin{1}) && length(varargin)==1
                obj.specs = [obj.specs, varargin{1}];
            elseif length(varargin)>1
                obj.specs = [obj.specs, varargin];
            else
                obj.specs = [obj.specs, {varargin{1}}];
            end

            %assign indices to the specifications
            for i = 1:numel(obj.specs)
                obj.specs{i}.id = i;
            end
        end

        function [rho, sperf] = perf_specs(obj, specs)
            %PERF_SPECS get the performance specifications and extract the
            %rho convergence rate 
            %
            % TODO: improve this part.
            %

            rho = 1;
            if nargin < 2
                sperf = obj.specs;
            else
                sperf = specs;
            end
            %extract the convergence rate
            %keep the stability spec if there is more than one
            %specification
            for i = 1:length(specs)
                if isa(specs{i}, 'spec_stability') 
                    if length(specs) > 1                        
                        sperf(i) = [];
                    end
                    rho = specs{i}.rho;
                end
            end

            for i =1:length(sperf)
                sperf{i}.id = i;
            end



        end

        function [sol] = solve_single(obj, arg, specs)
            %SOLVE_SINGLE Solve the program once
            %ADD_SPECIFICATIONS


            % warning('all', 'off')
            obj.specs = {};
            obj = obj.process_argument(arg);            
            obj = obj.add_specifications(specs);
            
            [vars, cons, objective, alg_psi, rho] = obj.build_program(); 

            % objective = obj.get_objective(vars);
            
            % objective = 0;
            [sol] = obj.run(vars, cons, objective);


            
            if sol.status == 0
                if obj.LMILAB
                    [vrec] = rec_vars(vars, sol.lmi_out);
                else
                    [vrec] = rec_vars(vars);
                end

                sol.vars = vrec;

                sol.rho = rho;
                sol.alg_psi = alg_psi;
                sol.objective = double(double(objective, sol.lmi_out));

                ncons = length(cons.lmim);

                if obj.config.recovery.blocks
                    sol.blocks = cell(ncons, 1);
                    for i = 1:ncons
                        sol.blocks{i} = double(double(cons.lmim(i), sol.lmi_out));
                    end
                end

                sol = obj.process_recovery(sol, sol.lmi_out, alg_psi);
                
            end

        end

        %% Bisection routines
        function [sol_best, vr] = bisect(obj, arg, specs, b_opts)
            %BISECT: perform bisection on a parameter. Minimization target
            %
            %sweep in options structure? 
            %
            %Input
            %
            % arg:      arguments for the routine (order for analysis, iqcs
            %           for synthesis)
            % specs:    performance specifications
            % (b_opts): bisection options (bisect_opts)
            %   param_range:    upper and lower bound for the parameter
            %   sweep_rho:      True:  sweep rho for the specification
            %                   False: sweep the bound for the specification
            %   spec_ind:       which specification in the list to sweep                   
            %
            %Output:
            %   sol_best: the best solution
            %   vr:       the range of the value at the optimal bisection


            sol_best = [];
            vr = [];
            %base relations

            %bisection options are separate from configuration options
            %should they be?
            if nargin < 4
                b_opts = bisect_opts;
            end
            
            

            %process the inputs and specifications
            obj = obj.process_argument(arg);
            
            obj.specs = [];
            obj = obj.add_specifications(specs);

            cons = obj.cons;
            specs = obj.specs;
            vars = obj.vars;

            % %form the plant
            % [iqc_data] = obj.iqc_op_all();
            % %TODO: modify this for different exponential weighting
            % [alg_psi, iqc_op, alg_loop]  = obj.sys.build_plant(iqc_data);
            % 
            % [vars, cons] = obj.lmi.create_vars(vars, cons, alg_psi, specs);            
            % 
            
            
            
            %take the initial step
            
            
            v = max(b_opts.val_range);
            vr = b_opts.val_range;
            WARM = false;


            
            %modify the specification use an oracle to do this
            
            % f = @(pcurr) bisect_inner(obj, pcurr, vars, cons, alg_psi, iqc_op, specs, b_opts);
            f = @(pcurr) bisect_inner(obj, pcurr, vars, cons, specs, b_opts);

            found_bound = [0, 0];
            sol0 = f(v);
            if sol0.status
                sol = sol0;
                vr = [vr(2), Inf];                
            else

                %TODO: implement warm starts

                % if WARM
                %     %search for upper bound
                % end

                while diff(vr) > b_opts.tol
                    %select midpoint (narrowing search)
                    v= sum(vr)/2;
                    [sol] = f(v);
                    valid_sol = (sol.status == 0);
                    if valid_sol;
                        %if successful, proceed with left interval 
                        vr(2) = v;
                        sol_best = sol;
                        v_best = v;
                    else
                        %if not successful, proceed with right interval
                        vr(1) = v;
                    end
                end
            end



            if ~isempty(sol_best) && sol_best.status == 0
                %ugly interface here
                vars_im = sol_best.vars;
                if obj.LMILAB
                    [vrec] = rec_vars(vars_im, sol.lmi_out);
                else
                    [vrec] = rec_vars(vars_im);
                end

                sol_best.vars = vrec;

                sol_best = obj.process_recovery(sol_best, sol_best.lmi_out, sol.alg_psi);
                sol_best.objective = double(double(sol_best.objective, sol_best.lmi_out));
            end

                        
        end

        function [spec_new] = modify_spec(obj, pcurr, spec_old, b_opts)
            %MODIFY_SPEC modify a specification in the bisection loop

            spec_new = spec_old;
            i = b_opts.spec_ind;
            if b_opts.bisect_rho && isa(spec_old{i}, 'spec_stability')
                spec_new{i}.rho = pcurr;
            else
                spec_new{i} = spec_new{i}.set_p(pcurr);                
            end


        end

        function obj = set_tol(obj, key, val)
            % SET_TOL set the tolerance of the LMI routines
            %
            % example: obj = obj.set_tol('G_max', 10)

            obj.lmi.tol = setfield(obj.lmi.tol, key, val);
        end

        function [sol] = bisect_inner(obj, pcurr, vars, cons, spec, b_opts)
            %BISECT_INNER: inner loop for bisection
            %run the program and process the solution
            
            % TODO: will need to redo constraint invocation due to the
            % exponential convergence implementation

            spec_curr = obj.modify_spec(pcurr, spec, b_opts);

            [vars, cons, objective, alg_psi, rho] = obj.build_program(spec_curr); 
            %form the plant
            % [iqc_data] = obj.iqc_op_all();
            %TODO: modify this for different exponential weighting
            % [alg_psi, iqc_op, alg_loop]  = obj.sys.build_plant(iqc_data);
            % 
            % [vars, cons] = obj.lmi.create_vars(vars, cons, alg_psi, spec_curr);            
            % 
            % [vars, cons, objective] = cons_dynamic(obj, vars, cons, alg_psi, iqc_op, spec_curr);
            
            
            [sol] = obj.run(vars, cons, objective);
            sol.objective = double(double(objective, sol.lmi_out));
            sol.rho = rho;
            sol.alg_psi = alg_psi;
            sol.vars = vars;
        end

        %% Constraint optimization

        %TODO: fmincon for the p2p objective

        %% Formation of Constraints and Plants
        function [iqc_data] = iqc_op_all(obj)
            %IQC_OP_ALL: all iqcs for the operators
            %
            %
            %useful for the build_plant routines
            iqc = {};
            m_same = [];
            ind_same = [];
            same_count = 0;

            for i = 1:length(obj.iqc_op)
                if ~obj.sys.op{i}.same
                    %block diagonal of the iqc
                    if isempty(iqc)
                        iqc = obj.iqc_op{i};
                    else
                        iqc = blkdiag(iqc, obj.iqc_op{i});
                    end
                    same_count = same_count + obj.iqc_op{i}.nw;

                else
                    %treat the m=L case separately
                    m_same = blkdiag(m_same, obj.iqc_op{i});

                    ind_same = [ind_same, same_count + (1:length(m_same))];
                    same_count = same_count + length(m_same);
                end            
            end

            iqc_data =struct('iqc', iqc, 'm_same', m_same, 'ind_same', ind_same, 'task', obj.task);

        end

        

        function [vars, cons, objective] = cons_dynamic(obj, vars, cons, alg_psi, iqc_op, specs)
            %CONS_DYNAMIC: form the dynamical dissipation relations for the
            %system (at the current set of specifications)
            %
            %Output:
            %   vars:       variables
            %   cons:       accumulated constraints
            %   objective:  single value to be minimized in inner loop (not
            %               the outer loop of bisection)

            [diss] = obj.index_specs(alg_psi, iqc_op, specs);
            ndiss = length(diss);

            %dissipation relations
            objective = 0;
            for i = 1:ndiss
                [cons, objective_curr] = obj.lmi.cons_dynamic(vars, cons, diss{i});                                 

                objective = objective + objective_curr;

                % sm = ssize(con_M, 1);
                % cons = append_lmi(cons, con_M - eye(sm)*obj.tol.M, obj.LMILAB);
                % cons = append_lmi(cons, con_M , obj.LMILAB);
            end

            %terminal cost/sign constraints
            % if obj.opts.impose_X
            %     %for infinite-horizon performance measures (l2 norm, h2
            %     %norm), terminal costs and sign constraints on the
            %     %storage function are not required. For finite-horizon
            %     %specifications (e.g. Invariance, peak-to-peak), they
            %     %are needed.                
            % 
            % end
        end

        %% getters and setters
        % function obj = set.config(obj, )
        % end


        
    end

    methods (Abstract)   
        process_argument(obj, arg); %inputs to run the routine
        process_recovery(obj, sol, lmi_out, alg_psi); %get the solution
        index_specs(obj, alg_psi, iqc_op, specs); %index the performance specifications
    end
end

