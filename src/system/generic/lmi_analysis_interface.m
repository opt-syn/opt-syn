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
       regcl;  %closed-loop regulator equation       
    end
    
    methods
        function obj = lmi_analysis_interface(sys, config)
            %LMI_ANALYSIS_INTERFACE Constructor
            %
            %Args:
            %   sys: algorithmic system
            %   config: configuration options
            %
            %Warning:
            %   Error: Regulator equation of system must be solvable,
            %          otherwise the algorithm is nonconvergent.

            obj@lmi_dispatch_interface(sys, config);

            %if the regulator equation fails, then the algorithm may
            %converge to a nonoptimal point.
            try
                obj.regcl = obj.reg.check_regulator();
            catch
                error('Failure of Regulator Equation: Convergence cannot be assured.')
            end
            
        end


        %% variable creation

        function [vars, cons] = create_vars(obj, vars, cons, alg_psi, specs)
            %CREATE_VARS create the variables for the problem
            %Args:            
            %   vars:  problem variables
            %   cons:  accumulated constraints
            %   alg_psi: generalized plant in analysis
            %   specs:  specifications
            %
            %Returns:
            %   vars:  problem variables
            %   cons:  accumulated constraints

            [vars.diss, cons] = obj.create_vars_storage(cons, alg_psi);
            [vars.spec, cons] = obj.create_vars_spec(cons, specs);
        end


        function [G, cons] = define_storage_G(obj, cons, alg_psi,  name)
            %DEFINE_STORAGE_G storage function for a specific subsystem
            %Args:            
            %   cons:  accumulated constraints
            %   alg_psi: generalized plant in analysis
            %   name:   name of the subsystem/object
            %
            %Returns:
            %   G:   storage function
            %   cons:  accumulated constraints

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
            %
            %Args:            
            %   cons:  accumulated constraints
            %   specs: cell of specification
            %
            %Returns:
            %   vars:  variables
            %   cons:  accumulated constraints
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
            %Args:
            %   G:  storage matrix
            %   cons:  accumulated constraints
            %   iqc_op: IQCs for the operators
            %
            %Returns:
            %   cons: accumulated constraints
            %   con_X: terminal constraint expression

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

    end

    methods (Abstract)
        %variable creation routines                              
    end
end

