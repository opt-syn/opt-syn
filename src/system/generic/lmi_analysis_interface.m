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
    
    % properties
        
        
    % end
    
    methods
        function obj = lmi_analysis_interface(sys, config)
            %LMI_ANALYSIS_INTERFACE Construct an instance of this class
            %   Detailed explanation goes here
            obj@lmi_dispatch_interface(sys, config);
        end


        %% variable creation

        function [vars, cons] = create_vars(obj, vars, cons, alg_psi, specs)
            %CREATE_VARS create the variables for the problem

            [vars.diss, cons] = obj.create_vars_storage(cons, alg_psi);
            [vars.spec, cons] = obj.create_vars_spec(cons, specs);
        end


        function [G, cons] = define_storage_G(obj, cons, alg_psi,  name)
            %DEFINE_STORAGE_G storage function for a specific subsystem
            
            if iscell(alg_psi)
                n = ssize(alg_psi{1}.A, 1); 
            else
                n = ssize(alg_psi.A, 1); 
            end
            G = lmim(['G', name], n, n, 'sym');

            

            if obj.config.tol.G_max < Inf    
                %issue in the bounding?
                cons = append_lmi(cons, obj.config.tol.G_max*eye(n)  - G, obj.config.LMILAB);
                cons = append_lmi(cons, obj.config.tol.G_max*eye(n)  + G, obj.config.LMILAB);                
            end            
        end

        function [vars_spec, cons] = create_vars_spec(obj, cons, specs)
            %CREATE_VARS_SPEC declare variables for the specifications

            %maybe put this somewhere else?
            %
            %right now the variables are in the (spec) object.
            nspec = length(specs);
            vars_spec = cell(nspec, 1);
            for i = 1:nspec
                [vars_spec{i}, cons] = specs{i}.create_vars(cons, [], obj.config);
            end
        end

        %% terminal constraints
        function [cons, con_X] = con_terminal(obj, G, cons, iqc_op)
            %CON_TERMINAL
            %terminal cost constraint (nonnegativity for the storage function G)
            %coupled positivity if the IQC has a terminal cost
            %
            %
            %
            X = iqc_op.X;
            

            nf = ssize(X);
            n = ssize(G, 1);
            Ef = [eye(nf); zeros(n-nf, nf)];

            X_f = Ef * X * Ef';
            con_X = G + X_f;

            sx = ssize(con_X, 1);

            if obj.config.gen.impose_X == 1
                cons = append_lmi(cons, con_X - eye(sx)*obj.config.tol.X, obj.config.LMILAB);
            elseif obj.config.gen.impose_X == 2
                cons = append_lmi(cons, G - eye(sx)*obj.config.tol.X, obj.config.LMILAB);
            end

        end
        
        %% common specification calls

        
        %function [cons, objective, con_M] = quad(obj, vars, cons, diss)
        %Quadratic performance (defined on a per-system basis)

        function [cons, objective, con_M] = stability(obj, vars, cons, diss)
            %STABILITY certification of exponential stability
            %
            %the supply function in the specification is empty,
            %so just call quadratic performance.


            [cons, objective, con_M] = obj.quad(vars, cons, diss);

        end

        function [cons, objective, con_M] = e2e(obj, vars, cons, diss)
            %E2E: energy to energy gain

            if diss.spec.target
                [cons, objective, con_M] = obj.e2e_target(vars, cons, diss);
            else
                %is a special case of quadratic performance
                [cons, objective, con_M] = obj.quad(vars, cons, diss);
            end           
        end

    end

    methods (Abstract)
        %variable creation routines        
        quad(obj, vars, cons, diss)               
    end
end

