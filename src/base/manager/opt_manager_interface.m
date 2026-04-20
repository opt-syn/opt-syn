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
        tol = struct('G', 1e-4, 'M', 1e-7, 'X', 1e-7);
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
    end
end

