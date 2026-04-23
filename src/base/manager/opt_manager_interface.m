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
                v_best = [];
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

