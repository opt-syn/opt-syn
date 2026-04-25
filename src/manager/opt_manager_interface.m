classdef (Abstract) opt_manager_interface
    %OPT_MANAGER_INTERFACE interface for the analysis and synthesis of
    %optimization/fixed point algorithms
    
    properties
        sys;  %system (opt_system type)
        cons = [];
        vars = {}; %specifications
        iqc_op = {};
        specs = {};        

        %numerics
        LMILAB = true;
        tol = struct('G_max', 100, ...
            'M', 1e-7, ...
            'X', 1e-7, ...
            'finite_l2', 1e3,...
            'dia', 1e-11);

        dispatch = [];
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

        function dis = select_dispatch(obj, sys)
            %SELECT_DISPATCH select the lmi dispatch routines
            %
            % dispatch: handlers for the LMIs
            %
            %classify the type of the system

            %TODO: need to implement this
            if isa(sys, 'opt_system')
                %TODO: working on it
                dis = lmi_dispatch_lti(sys);
            elseif isa(sys, 'opt_system_periodic')
                %TODO: not yet implemented
                dis = lmi_dispatch_periodic(sys);
            elseif isa(sys, 'opt_system_switched')
                
                if all((sys.adj==0) + (sys.adj==1), 'all')
                    %robust switching
                    %TODO: not yet implemented
                    dis = lmi_dispatch_switched(sys);
                else
                    %TODO: not yet implemented
                    %stochastic: markof jump linear system
                    dis = lmi_dispatch_mjls(sys);
                end                
            elseif isa(sys, 'opt_system_lpv_poly')
                %maybe? or go to switched system
                dis = lmi_dispatch_lpv_poly(sys);
            elseif isa(sys, 'opt_system_lpv')
                dis = lmi_dispatch_lpv_lfr(sys);
                %TODO: not yet implemented
            end

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

            if STATUS == 0
                if obj.LMILAB
                    [vrec] = rec_vars(vars, lmi_out);
                else
                    [vrec] = rec_vars(vars);
                end

                sol.vars = vrec;

                sol = obj.process_recovery(sol, lmi_out);
            end
            
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
            if iscell(varargin{1})
                obj.specs = [obj.specs, varargin{1}];
            else
                obj.specs = [obj.specs, varargin];
            end
        end

        function [vars_spec, cons] = create_vars_spec(obj, specs, cons)
            %CREATE_VARS_SPEC declare variables for the specifications
            nspec = length(specs);
            vars_spec = cell(nspec, 1);
            for i = 1:nspec
                [vars_spec{i}, cons] = specs{i}.create_vars(cons);
            end

        end

        function [sol] = solve_single(obj, order, specs)
            %SOLVE_SINGLE Solve the program once

            if ~iscell(order)
                order = {order};               
            end

            % warning('all', 'off')
            obj = obj.oracle_order(order);
            obj = obj.add_specifications(specs);

            [vars, cons] = obj.build_program(); 

            objective = obj.get_objective(vars);
            
            [sol] = obj.run(vars, cons, objective);

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

            vars = obj.vars;
            [alg_psi, iqc_op, alg_loop]  = obj.build_plant(cons);
            [vars.diss, cons] = obj.create_vars_storage(alg_psi, cons);            
            [vars.spec, cons] = obj.create_vars_spec(specs, cons);
            
            
            
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
                    sol = f(v);
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

        function [sol] = bisect_inner(obj, pcurr, vars, cons, alg_psi, iqc_op, spec, b_opts)
            %BISECT_INNER: inner loop for bisection
            %run the program and process the solution
            objective = obj.get_objective(vars);

            spec_curr = obj.modify_spec(pcurr, spec, b_opts);

            [vars, cons] = build_dissipation(obj, vars, cons, alg_psi, iqc_op, spec_curr);
            [sol] = obj.run(vars, cons, objective);
            
        end

        %% Indexing routines for systems
                function [diss] = index_specs(obj, alg_psi, iqc_op, specs)

            %INDEX_SPECS:  index into the performance specifications
            %
            %
            %now index alg_psi into its performance specifications
            if nargin < 4
                specs = obj.specs;
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


                %need to permute the entries of Mdiag for the partition
                n1 = iqc_op.np;
                m1 = iqc_op.nq;
                n2 = length(sp.iwp);
                m2 = length(sp.izp);
                [M] = outer_blkdiag(iqc_op.M, sp.supply, n1, m1, n2, m2);
                % Mdiag = blkdiag(iqc_op.M, sp.supply);




                %TODO: this may run into trouble if one entry has an X.
                %performance with dynamic multipliers?
            
                diss{i} = struct('plant', alg_screen, 'M', M, 'X', iqc_op.X, ...
                    'spec', sp);
            end

        end


        
    end

    methods (Abstract)
        oracle_order(obj)
        build_dissipation(obj)
        build_program(obj)
        build_plant(obj, cons)
        get_objective(obj)
        process_recovery(obj, sol);
    end
end

