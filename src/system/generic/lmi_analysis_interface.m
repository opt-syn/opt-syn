classdef lmi_analysis_interface < lmi_dispatch_interface
    %LMI_ANALYSIS_INTERFACE 
    %Linear Matrix Inequality constraints for analysis of algorithmic
    %interconnections.
    %
    %
    %this is overridden by specialized analysis routines for system types:
    %   lti
    %   periodic
    %   switched robust
    %   switched jump
    
    properties
        tol = struct('G_max', 100 ... %upper bound on norm of storage matrix
            );
        LMILAB = 1;
    end
    
    methods
        function obj = lmi_analysis_interface(sys)
            %LMI_ANALYSIS_INTERFACE Construct an instance of this class
            %   Detailed explanation goes here
            obj@lmi_dispatch_interface(sys);
        end


        %% variable creation

        function [vars, cons] = create_vars(obj, vars, cons, alg_psi, specs)
            %CREATE_VARS create the variables for the problem

            [vars.diss, cons] = obj.create_vars_storage(alg_psi, cons);
            [vars.spec, cons] = obj.create_vars_spec(specs, cons);
        end

        function [vars_diss, cons]= create_vars_storage(obj, alg_psi, cons, name)
            %create_vars_storage create variables for the dissipation
            %constraints

            if nargin < 4
                name = [];
            end

            n = length(alg_psi.A);
            G = lmim(['G', name], n, n, 'sym');

            vars_diss= struct('G', G);
 
            if obj.tol.G_max < Inf    
                %issue in the bounding?
                cons = append_lmi(cons, obj.tol.G_max*eye(n)  - G, obj.LMILAB);
                cons = append_lmi(cons, obj.tol.G_max*eye(n)  + G, obj.LMILAB);

                %lmim complains that the dimensions are wrong here.
                % cons = append_lmi(cons, ga*eye(n)  - G, obj.LMILAB);
                % cons = append_lmi(cons, ga*eye(n)  + G, obj.LMILAB);
            end
        end

        function [vars_spec, cons] = create_vars_spec(obj, specs, cons)
            %CREATE_VARS_SPEC declare variables for the specifications

            %maybe put this somewhere else?
            %
            %right now the variables are in the (spec) object.
            nspec = length(specs);
            vars_spec = cell(nspec, 1);
            for i = 1:nspec
                [vars_spec{i}, cons] = specs{i}.create_vars(cons);
            end

        end

    end

    % methods (Abstract)
    %     %variable creation routines        
    %     create_vars_spec(obj, specs, cons)
    % end
end

