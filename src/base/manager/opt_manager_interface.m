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
            'finite_l2', 1e3);

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

        function [sol] = run(obj,  vars, cons, objective)

            sol = struct;
            %run the program
            if obj.LMILAB
                if ~isnumeric(objective)
                    lmis(cons, objective, 'c');
                end
                [lmi_out,info_out]=lmisolve(cons);
                
                STATUS = lmi_out.status;
                sol.dia = lmi_out.dia;
                sol.info = info_out;
                ncons = length(cons.lmim);
                sol.blocks = cell(ncons, 1);
                sol.eb = zeros(ncons, 1);
                for i = 1:length(cons.lmim)
                    sol.blocks{i} = -double(double(cons.lmim(i), lmi_out));
                    sol.eb(i) = min(eig(sol.blocks{i}));
                end
                
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
            obj.specs = [obj.specs, varargin];
        end
        
    end

    methods (Abstract)
        build_dissipation(obj)
        build_program(obj)
        build_plant(obj)
        process_recovery(obj, sol);
    end
end

