classdef (Abstract) opt_manager_interface
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
        
        
        %numerics
        LMILAB = true;
        tol = struct('dia', 1e-11);          %tolerance for acceptable solution 

        lmi = [];
        %other options
        opts = struct('impose_X', true)
    end
    
    methods
        function obj = opt_manager_interface(sys)
            %OPT_MANAGER_INTERFACE
            obj.sys = sys;
            nop = length(obj.sys.op);
            obj.iqc_op = cell(nop, 1);
            obj.vars = struct('op', []);
            obj.vars.op =cell(nop, 1);            

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

            lmi_handler = lmi_hand(sys);

        end



        function [vars, cons, objective] = build_program(obj, specs)
            %BUILD_PROGRAM set up the algorithm analysis problem

            if nargin < 2
                specs = obj.specs;
            end

            cons = obj.cons;
            vars = obj.vars;

            % [diss, cons] = obj.build_plant(cons);
            %alg_loop: used for debugging. The algorithm after signal 
            % transformationsbefore, but before cascade by the filters    

            [iqc_data] = obj.iqc_op_all();
            [alg_psi, iqc_op, alg_loop] = obj.sys.build_plant(iqc_data);

            [vars, cons] = obj.lmi.create_vars(vars, cons, alg_psi, specs);

            %the dissipation can change
            [vars, cons, objective] = obj.cons_dynamic(vars, cons, alg_psi, iqc_op, specs);

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
                STATUS = (lmi_out.status || (sol.dia > obj.tol.dia));
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
                opt = sdpsettings('verbose', p.opts.verbose, 'solver', p.opts.solver);
            
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


        function [sol] = solve_single(obj, order, specs)
            %SOLVE_SINGLE Solve the program once

            %ADD_SPECIFICATIONS
            if ~iscell(order)
                order0 =order;
                nop = length(obj.sys.op);
                order = cell(nop, 1);
                for i = 1:nop
                    order{i}  = order0;
                end
            end

            % warning('all', 'off')
            obj = obj.oracle_order(order);
            obj = obj.add_specifications(specs);

            [vars, cons, objective] = obj.build_program(); 

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

                sol = obj.process_recovery(sol, sol.lmi_out);
                sol.objective = double(double(objective, sol.lmi_out));
            end

        end

        %% Bisection routines
        function [sol_best, vr] = bisect(obj, order, specs, b_opts)
            %BISECT: perform bisection on a parameter. Minimization target
            %
            %sweep in options structure? 
            %
            %Input
            %
            % (b_opts): bisection options (bisect_opts)
            %   param_range:    upper and lower bound for the parameter
            %   sweep_rho:      True:  sweep rho for the specification
            %                   False: sweep the bound for the specification
            %   spec_ind:       which specification in the list to sweep                   
            %
            %Output:
            %   sol_best: the best solution
            %   vr:       the range of the value at the optimal bisection

            if ~iscell(order)
                order = {order};               
            end
            if nargin < 4
                b_opts = bisect_opts;
            end

            sol_best = [];
            vr = [];
            %base relations
            
            obj = obj.oracle_order(order);
            obj = obj.add_specifications(specs);
            cons = obj.cons;
            specs = obj.specs;

            vars = obj.vars;
            [iqc_data] = obj.iqc_op_all();
            [alg_psi, iqc_op, alg_loop]  = obj.sys.build_plant(iqc_data);
            [vars, cons] = obj.lmi.create_vars(vars, cons, alg_psi, specs);            
            
            
            
            
            %take the initial step
            
            % if isempty(b_opts.val_init)
                v = max(b_opts.val_range);
                vr = b_opts.val_range;
                WARM = false;
            % else
            %     v = b_opts.val_init;
            %     vr = [1 - b_opts.warm_factor, 1 + b_opts.warm_factor] * v;
            %     WARM = true;
            % end


            
            %modify the specification use an oracle to do this
            
            f = @(pcurr) bisect_inner(obj, pcurr, vars, cons, alg_psi, iqc_op, specs, b_opts);

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

            if sol.status == 0
                if obj.LMILAB
                    [vrec] = rec_vars(vars, sol.lmi_out);
                else
                    [vrec] = rec_vars(vars);
                end

                sol.vars = vrec;

                sol = obj.process_recovery(sol, sol.lmi_out);
                sol.objective = double(double(sol.objective, sol.lmi_out));
            end

                        
        end

        function [spec_new] = modify_spec(obj, pcurr, spec_old, b_opts)
            %MODIFY_SPEC modify a specification in the bisection loop

            spec_new = spec_old;
            i = b_opts.spec_ind;
            if b_opts.bisect_rho
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

        function [sol] = bisect_inner(obj, pcurr, vars, cons, alg_psi, iqc_op, spec, b_opts)
            %BISECT_INNER: inner loop for bisection
            %run the program and process the solution
            

            spec_curr = obj.modify_spec(pcurr, spec, b_opts);

            [vars, cons, objective] = cons_dynamic(obj, vars, cons, alg_psi, iqc_op, spec_curr);
            [sol] = obj.run(vars, cons, objective);
            sol.objective = double(double(objective, sol.lmi_out));
            
            
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

        function [diss] = index_specs(obj, alg_psi, iqc_op, specs)

            %INDEX_SPECS:  index into the performance specifications
            %
            %
            %   diss:   structure describing the problem
            %       plant:  system to control
            %       spec:   performance specification           
            %       target: whether the performance measure should be optimized
            %               true:  soft constraint (e.g. Schur complement
            %                                       formulation)
            %               false: hard constraint
            
            %TODO: figure out synthesis here
            
            if nargin < 4
                specs = obj.specs;
            end

            if nargin < 5
                target_ind = 0;
            end

            diss = cell(length(specs), 1);
            %determine the indices for each performance specification
            for i = 1:length(specs)
                
                      
                sp = specs{i};
                iwp_iqc = (1:(iqc_op.nw))';
                ir_iqc_first = (1:(iqc_op.np))';


                count_iqc_in = (iqc_op.nw);
                count_iqc_out = (iqc_op.np);

                if isempty(sp.izp) || isempty(sp.iwp)
                    ir_iqc_first_r =[];
                    iw_iqc_first_r = [];
                else
                    iw_iqc_first_r = count_iqc_in + (1:sp.iwp);
                    count_iqc_in = count_iqc_in + sp.iwp;

                    ir_iqc_first_r = count_iqc_out  + (1:sp.izp);
                    count_iqc_out = count_iqc_out + sp.izp;
                end

                iwp_iqc = [iwp_iqc; iw_iqc_first_r];

                ir_iqc = [ir_iqc_first; ir_iqc_first_r];
                ir_iqc = [ir_iqc; ir_iqc + (iqc_op.np + obj.sys.P.nwp )];

                sp_ind_w = iwp_iqc;
                sp_ind_r = ir_iqc;
                % nww = length(sp_ind_w);
                % nwr = length(sp_ind_r);

                %TODO: allow for synthesis here
                [nwr, nww] = ssize(alg_psi.D);
                %enforce squareness in the performance specs?
                E_w = full(sparse(1:length(sp_ind_w), sp_ind_w, ones(1, length(sp_ind_w)), length(sp_ind_w), nww));
                E_r = full(sparse(1:length(sp_ind_r), sp_ind_r, ones(1, length(sp_ind_r)), length(sp_ind_r), nwr));


                %nonminimal representation
                alg_screen = E_r * alg_psi * E_w;


                diss{i} = struct('plant', alg_screen, 'iqc_rob', iqc_op, ...
                    'spec', sp);
                % %need to permute the entries of Mdiag for the partition





                %TODO: this may run into trouble if one entry has an X.
                %performance with dynamic multipliers?
            
                % diss{i} = struct('plant', alg_screen, 'M', M, 'X', iqc_op.X, ...
                    % 'spec', sp);
            end

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


        
    end

    methods (Abstract)        
        oracle_order(obj)                        
        process_recovery(obj, sol);
    end
end

