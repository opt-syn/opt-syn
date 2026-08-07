classdef (Abstract) opt_manager_interface < handle
    %OPT_MANAGER_INTERFACE interface for the analysis and synthesis of
    %optimization/inclusion algorithms

    
    %inheritance:
    % opt_analysis  < opt_manager_interface 
    % opt_synthesis < opt_manager_interface 
    
    properties
        sys;  %system (opt_system type)        
        cons = []; %accumualted constraints
        vars = {}; %variables of the problem.
        iqc_op = {}; %iqcs for the operators
        specs = {};  %performance specifications
        
        config = []; %configuration options (opt_config)

        lmi = []; % (lmi_dispatch) the lmi handler 
        %other options
        task = 'generic'; %analysis or synthesis? This is an override
        sys_orig; %(opt_system) original system, before any performance-based modifications 
    end

    
    methods
        function obj = opt_manager_interface(sys, config)
            %OPT_MANAGER_INTERFACE Constructor
            % Args:
            %   sys: algorithmic system
            %   config: configuration options
            
            if nargin < 2
                config = opt_config();
            end

            obj.sys = sys;
            obj.sys_orig = sys;
            nop = length(obj.sys.op);
            obj.iqc_op = cell(nop, 1);
            obj.vars = struct('op', []);
            obj.vars.op =cell(nop, 1);   

            obj.config = config;
            

        end

        function lmi_handler = select_lmi(obj, sys)
            %SELECT_LMI select the lmi routines based on the system type
            %
            %Args:
            %   sys:    type of system (e.g. opt_system_switched)
            %   routine:    'analysis' or 'synthesis'
            %Returns:
            %   lmi_handler: the lmi object for the specific dynamics            

            tp = sys.get_type();
            clname = ['lmi_', obj.task, '_', tp];

            
            %hard-code override for reduced-order synthesis
            if strcmp(clname, 'lmi_synthesis_lti') && obj.config.syn.reduced_order
                clname = 'lmi_synthesis_lti_reduced_order';
            end

            if strcmp(clname, 'lmi_synthesis_periodic_orbit') && obj.config.syn.reduced_order
                clname = 'lmi_synthesis_periodic_orbit_reduced_order';
            end

            lmi_hand = str2func(clname);
            

            lmi_handler = lmi_hand(sys, obj.config);

        end



        function [vars, cons, objective, alg_psi, diss] = build_program(obj, specs)
            %BUILD_PROGRAM set up the algorithm analysis or synthesis problem
            %Args:
            %   specs: specifications            
            %Returns:
            %   vars:   variables of the problem        
            %   cons (lmibl):   accumulated constraints
            %   objective: objective to minimize
            %   alg_psi (genplant/genplantpoly):    generalized plant,
            %   before internal model
            %   diss:   current dissipation inequality
                      

            if nargin < 2
                specs = obj.specs;
            end

            cons = obj.cons;
            vars = obj.vars;
            
            %alg_loop: used for debugging. The algorithm after signal 
            % transformationsbefore, but before cascade by the filters    

            [sperf, ERGODIC] = obj.perf_specs(specs);

            
            %override for same rho
            if obj.config.gen.same_rho

                %compensate: the current code assumes that the IQC is
                %rho-weighted. 
                common_rho = obj.get_common_rho(sperf);                

                %enforce all attributes to have the same rho
                ind_stab = 0;
                for i = 1:length(sperf)
                    sperf{i}.rho = common_rho;
                    if isa(sperf{i}, 'spec_stability')
                        ind_stab = i;
                    end
                end

                %drop the exponential stability specification if it is
                %alone                
                if (ind_stab > 0) && (length(sperf)>1)
                    sperf(ind_stab) = [];
                    for i = ind_stab:length(sperf)
                        sperf{i}.id = sperf{i}.id - 1;
                    end
                end



                if strcmp(obj.task, 'synthesis')
                    % weight each IQC from analysis, then factor each IQC

                    %this needs to be done individually for each rho
                    iqc_op_factored = obj.iqc_op_ana;
                    for i = 1:length(iqc_op_factored)
                        if isnumeric(obj.iqc_op_ana{i})
                            iqc_op_factored{i} = obj.iqc_op_ana{i};
                        else
                            iqccurr = rhotrafo(obj.iqc_op_ana{i}, common_rho);
                        
                            iqc_op_factored{i} = iqccurr.factor();
                        end
                    end
                    obj.iqc_op = iqc_op_factored;
                    [iqc_data] = obj.iqc_op_all();
                    iqc_data.iqc = iqc_data.iqc.rhotrafo(1/common_rho);
                else
                    [iqc_data] = obj.iqc_op_all();
                end
                
                
            else
                %either analysis or the nice situation
                [iqc_data] = obj.iqc_op_all();
            end
            iqc_data.ERGODIC = ERGODIC;

            
            [alg_psi, iqc_op, alg_loop] = obj.sys.build_plant(iqc_data);



            iqc_data.iqc_op = iqc_op;
            [vars, cons] = obj.lmi.create_vars(vars, cons, alg_psi, sperf);

            %the dissipation can change
            [vars, cons, objective, diss] = obj.cons_dynamic(vars, cons, alg_psi, iqc_data, sperf);

            %add spreading constraint
            cons = obj.lmi.con_spread(cons, vars);

        end

        function verdict = LMILAB(obj)
            %is LMILAB used?
            %Returns:
            %   verdict (bool): 
            verdict = obj.config.LMILAB();
        end

        %% acquire solutions
        function [sol] = run(obj,  vars, cons, objective)
            %RUN: run the program
            % Args:
            %   vars:   variables of the problem        
            %   cons:       accumulated constraints
            %   objective:  target to minimize
            % Returns
            %   sol: solution structure

            
            sol = opt_solution(obj.task);
            
            
            if obj.LMILAB
                if isnumeric(objective)
                    OBJECTIVE = false;                                     
                else
                    cons = lmis(cons, objective, 'c');
                    OBJECTIVE = true;
                    
                end
                
                
                [lmi_out,info_out]=lmisolve(cons);
                
                

                STATUS = (lmi_out.status || (lmi_out.dia(1+OBJECTIVE) > obj.config.tol.dia));
                sol.info = info_out;
                sol.info.cons = cons;

                
                if OBJECTIVE
                    if STATUS
                        sol.dia = lmi_out.dia;
                        sol.objective = Inf;
                    else
                        sol.objective = lmi_out.dia(1);
                        sol.dia = lmi_out.dia(2);
                    end
                    
                else
                    sol.dia = lmi_out.dia;
                    sol.objective = 0;
                end

                %explicitly forming the blocks in recovery takes a long
                %time. Why is this the case?

                
                if obj.config.recovery.blocks
                    ncons = length(cons.lmim);
                    sol.recovery = struct;
                    sol.recovery.blocks = cell(ncons, 1);
                    sol.recovery.eb = zeros(ncons, 1);
                    for i = 1:length(cons.lmim)
                        sol.recovery.blocks{i} = -double(double(cons.lmim(i), lmi_out));
                        sol.recovery.eb(i) = min(eig(sol.recovery.blocks{i}));
                    end
                end
                
            else %YALMIP
                %warning: not implemented
                opt = sdpsettings('verbose', obj.gen.config.verbose, 'solver', obj.config.gen.solver);
            
                t = optimize(lmi, objective, opt);
                
                STATUS = t.problem;
            end
            sol.status = STATUS;
            sol.info.lmi_out = lmi_out;          
        end


        function obj = purge(obj);
            %PURGE: get rid of problem description
            %Returns:
            %   obj: cleaned up manager
            obj.cons = [];
            obj.vars = {};
            obj.iqc = {};
            obj.specs = {};
        end
        
        function obj = scan_specifications(obj, varargin)
            %concatenate the new performance specifications
            %Args:
            %   varargin: new specifications to add
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

                %add extra performance channels for warranted
                %specifications (e.g. ergodic)
                if isa(obj.specs{i}, 'spec_ergodic')
                    %requires duality gap condition, add it here.

                    
                    %propagate through the updates
                    [obj.specs{i}, obj.sys] = obj.specs{i}.augment_ergodic(obj.sys);
                    obj.lmi.sys = obj.sys;
                    obj.lmi.reg.sys = obj.sys;
                end

            end
        end


        function rho = get_common_rho(obj, specs)
            %GET_COMMON_RHO get the common rho in the case of the same rho in all
            %performance specifications. required to use noncausal
            %multipliers
            %
            %Args:
            %   specs (cell):  array of specifications
            %
            %Returns:
            %   rho (double): the common rho 
            rho = 1;

            if obj.config.gen.same_rho
                for i = 1:length(specs)
                    if isa(specs{i}, 'spec_stability')
                        rho = specs{i}.rho;
                    end
                end
            end
        end

        function [sperf, ERGODIC] = perf_specs(obj, specs)
            %PERF_SPECS index the  performance specifications 
            %
            %Args:
            %   specs (cell):  performance specifications
            %
            %Returns:
            %   sperf: the specification cell 
            %   ERGODIC: is ergodic convergence required


            if nargin < 2
                sperf = obj.specs;
            else
                sperf = specs;
            end
            %extract the convergence rate
            %keep the stability spec if there is more than one
            %specification
            ERGODIC = false;
            for i = 1:length(specs)
                if isa(specs{i}, 'spec_ergodic')
                    ERGODIC = true;
                end
            end

            for i =1:length(sperf)
                sperf{i}.id = i;                
            end



        end

        function [sol] = solve_single(obj, arg, specs)
            %SOLVE_SINGLE Solve the program once
            %
            %Args:
            %   arg:    order (analysis) or iqc (synthesis)
            %   specs:  specification cell


            % warning('all', 'off')
            if nargin < 3 || isempty(specs)
                specs = {spec_stability(1)};
            end
            if ~iscell(specs)
                specs = {specs};
            end
            obj.specs = {};
            obj = obj.process_argument(arg);            
            obj = obj.scan_specifications(specs);
            
            [vars, cons, objective, alg_psi, diss] = obj.build_program(); 

            % objective = obj.get_objective(vars);
            

            % objective = 0;
            [sol] = obj.run(vars, cons, objective);


            
            if sol.status == 0
                if obj.LMILAB
                    [vrec] = rec_vars(vars, sol.info.lmi_out);
                else
                    [vrec] = rec_vars(vars);
                end

                sol.vars = vrec;

                for i = 1:length(obj.specs)
                    if isa(obj.specs{i}, 'spec_stability')
                        sol.rho = obj.specs{i}.rho;
                    end
                end
                sol.cert.alg_psi = alg_psi;
                sol.objective = double(double(objective, sol.info.lmi_out));

                ncons = length(cons.lmim);

                sol = obj.process_recovery(sol, sol.info.lmi_out, alg_psi, diss);
                
            end

        end

        %% Bisection routines
        function [sol_best, vr] = bisect(obj, arg, specs)
            %BISECT: perform bisection on a parameter to minimize an objective.
            %
            % Args:
            %   arg:     arguments for the routine (order for analysis, iqcs for synthesis)
            %   specs:    (cell) performance specifications            
            % Returns:
            %   sol_best: the best solution
            %   vr:       the range of the value at the optimal bisection


            sol_best = [];
 
            %base relations

            %bisection options are separate from configuration options
            %should they be?

            if nargin < 2
                arg = [];
            end
            if nargin < 3 || isempty(specs)
                specs = {spec_stability(1)};
            end

            
            %process the inputs and specifications
            obj = obj.process_argument(arg);
            
            obj.specs = [];
            obj = obj.scan_specifications(specs);

            cons = obj.cons;
            specs = obj.specs;
            vars = obj.vars;

 
            
            
            %take the initial step
            
            
            v = max(obj.config.bisect.val_range);
            vr = obj.config.bisect.val_range;
            
            
            %TODO: implement warm start logic
            WARM = false;


            
            %modify the specification use an oracle to do this
            
            % f = @(pcurr) bisect_inner(obj, pcurr, vars, cons, alg_psi, iqc_op, specs, b_opts);
            f = @(pcurr) bisect_inner(obj, pcurr, vars, cons, specs);

            % found_bound = [0, 0];
            sol0 = f(v);
            if sol0.status
                sol_best = sol0;
                vr = [vr(2), Inf];                
            else

                %TODO: implement warm starts

                % if WARM
                %     %search for upper bound
                % end

                while diff(vr) > obj.config.bisect.tol
                    %select midpoint (narrowing search)
                    v= sum(vr)/2;
                    [sol] = f(v);
                    valid_sol = (sol.status == 0);
                    if valid_sol;
                        %if successful, proceed with left interval 
                        vr(2) = v;
                        sol_best = sol;
                        % v_best = v;
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
                    [vrec] = rec_vars(vars_im, sol_best.info.lmi_out);
                else
                    [vrec] = rec_vars(vars_im);
                end

                sol_best.vars = vrec;

                sol_best = obj.process_recovery(sol_best, sol_best.info.lmi_out, sol_best.cert.alg_psi, sol_best.cert.diss);
                sol_best.objective = double(double(sol_best.objective, sol_best.info.lmi_out));
            end

                        
        end

        function [spec_new] = modify_spec(obj, pcurr, spec_old, b_opts)
            %MODIFY_SPEC modify a specification in the bisection loop
            %
            %Args:
            %   pcurr: current value of the parameter to set
            %   spec_old:   specification to update            
            %   
            %Returns:
            %   spec_new: updated specification
            spec_new = spec_old;
            i = obj.config.bisect.spec_ind;
            if obj.config.bisect.bisect_rho || isa(spec_old{i}, 'spec_stability')
                spec_new{i}.rho = pcurr;
            else
                spec_new{i} = spec_new{i}.set_p(pcurr);                
            end


        end

        function obj = set_tol(obj, key, val)
            % SET_TOL set the tolerance of the LMI routines
            %
            % Example: obj = obj.set_tol('G_max', 10)

            obj.lmi.tol = setfield(obj.lmi.tol, key, val);
        end

        function [sol] = bisect_inner(obj, pcurr, vars, cons, spec, b_opts)
            %BISECT_INNER: inner loop for bisection
            %run the program and process the solution
            %Args:
            %   pcurr: current value of the parameter to set
            %   vars:   variables of problem
            %   cons:  accumulated constraints
            %   spec:   specifications            
            %   
            %Returns:
            %   sol: solution structure

            % TODO: will need to redo constraint invocation due to the
            % exponential convergence implementation

            spec_curr = obj.modify_spec(pcurr, spec);

            [vars, cons, objective, alg_psi, diss] = obj.build_program(spec_curr); 
 
            [sol] = obj.run(vars, cons, objective);

            if obj.config.bisect.bisect_rho
                sol.rho = pcurr;
            else
                for i = 1:length(spec_curr)
                    if isa(spec_curr{i}, 'spec_stability')
                        sol.rho = spec_curr{i}.rho;
                    end
                end

            end
            
            
            sol.spec = spec_curr;
            sol.cert.alg_psi = alg_psi;
            sol.vars = vars;
            sol.cert.diss = diss;
        end

        %% Constraint optimization

        %TODO: fmincon for the p2p objective

        %% Formation of Constraints and Plants
        function [iqc_data] = iqc_op_all(obj, iqc_op)
            %IQC_OP_ALL: all iqcs for the operators
            %
            % Returns:
            %   iqc_data (iqc_data_container): information for the iqcs
            
            if nargin< 2
                iqc_op = obj.iqc_op;
            end
            
            %useful for the build_plant routines
            iqc = {};
            m_same = [];
            ind_same = [];
            same_count = 0;            

            for i = 1:length(iqc_op)
                c = obj.sys.op{i}.c;
                if ~obj.sys.op{i}.same
                    %block diagonal of the iqc
                    if isempty(iqc)
                        iqc = iqc_op{i};
                    else
                        iqc = blkdiag(iqc, iqc_op{i});
                    end
                    same_count = same_count + (c);

                else
                    %treat the m=L case separately
                    m_same = blkdiag(m_same, kron(iqc_op{i}, eye(c)));

                    ind_same = [ind_same, same_count + (1:(c))];
                    same_count = same_count + length(m_same);
                end            
            end

            
            iqc_data =iqc_data_container;
            iqc_data.iqc = iqc;
            iqc_data.m_same = m_same;
            iqc_data.ind_same = ind_same;
            iqc_data.task =  obj.task;    

            iqc_data.augmented = false;
            if strcmp(obj.sys.type, 'periodic_orbit')
                iqc_data.rotate = true;
            end

        end

        

        function [vars, cons, objective, diss] = cons_dynamic(obj, vars, cons, alg_psi, iqc_data, specs)
            %CONS_DYNAMIC: form the dynamical dissipation relations for the
            %system (at the current set of specifications)
            %
            %Returns:
            %   vars:   variables of the problem        
            %   cons:       accumulated constraints
            %   objective:  single value to be minimized in inner loop (not
            %               the outer loop of bisection)            

            [diss] = obj.index_specs(alg_psi, iqc_data, specs);
            ndiss = length(diss);

            %dissipation relations
            objective = 0;
            for i = 1:ndiss
                [vars, cons, objective_curr] = obj.lmi.cons_dynamic(vars, cons, diss{i});                                 

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

        
    end

    methods (Abstract)   
        process_argument(obj, arg); %inputs to run the routine
        process_recovery(obj, sol, lmi_out, alg_psi, diss); %get the solution
        index_specs(obj, alg_psi, iqc_op, specs); %index the performance specifications
    end
end

