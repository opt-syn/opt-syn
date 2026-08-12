classdef lmi_analysis_interface < lmi_dispatch_interface
    %LMI_ANALYSIS_INTERFACE 
    %Linear Matrix Inequality constraints for analysis of algorithmic
    %interconnections.  
    %
    %
    %this is overridden by specialized analysis routines for system types.   
    
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

            if isempty(sys.K) || (iscell(sys.K) && isempty(sys.K{1}))

            
                error('Analysis: No controller is provided in the System.')
                
            end

            %if the regulator equation fails, then the algorithm may
            %converge to a nonoptimal point.
            try
                obj.regcl = obj.reg.check_regulator();
            catch
                error('Analysis: Failure of Regulator Equation: Convergence cannot be assured.')
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

            if isempty(iqc_op)
                %no uncertainty is present
                sx = ssize(G, 1);
                cons = append_lmi(cons, G - eye(sx)*obj.config.tol.X, obj.config.LMILAB);
            else
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

        end

        function objective = default_objective(obj, vars)
            %  create an objective for mincx() if one doesn't already
            %  exist.
            %                       
            %Args:
            %   vars:  variables of the problem
            %Returns:                 
            %   objective:   the objective to minimize


            if iscell(vars.diss)
                objective = 0;
                for i = 1:length(vars.diss)
                    objective = objective + trace(vars.diss.G);
                end
            else
                objective = trace(vars.diss.G);
            end

        end

        %% recovery
        function sol = process_recovery(obj, sol, lmi_out, alg_psi, diss)
            %recover the controller
            %Args:
            %   sol: solution structure
            %   lmi_out: output from solver
            %   alg_psi:   the filtered algorithmic interconnection
            %   diss (diss_data):   structure describing the dissipation constraint            
            %
            %Returns:  
            %   sol: solution structure


            %system-dependent overrides
        end
        
        %% common specification calls
        function [plant_rob, plant_perf] = partition_perf(obj, diss)
            %PARTITION_PERF partition the robust and performance channels
            %
            %Args:
            %   diss (diss_data):   current dissipation constraint
            %
            %Returns:
            %   plant_rob (sdpss):   plant with robust outputs
            %   plant_perf (sdpss):  plant with performance outputs
            nz_all = diss.iqc_rob.nq + diss.iqc_rob.nq;
            nwp = length(diss.spec.iwp);
            nzp = length(diss.spec.izp);

            I = eye(nz_all + nwp + nzp);
            I_robust = I(1:nz_all, :);
            I_perf = I((nz_all+1):end, :);

            plant_rob = I_robust * diss.plant;
            
            if nwp + nzp > 0
                plant_perf = I_perf * diss.plant;
            else
                plant_perf = [];
            end
              

        end

        function [plant_rob, plant_perf] = partition_perf_zp(obj, diss)
            %PARTITION_PERF_zp partition the robust and performance channels
            %ignore the performance inputs wp
            %
            %Args:
            %   diss (diss_data):   current dissipation constraint
            %
            %Returns:
            %   plant_rob (sdpss):   plant with robust outputs
            %   plant_perf (sdpss):  plant with performance outputs
            
            %input indexer
            nwp = diss.spec.nwp;
            nw = ssize(diss.plant.D, 2) - nwp;
            E_rob = screen_system([nw, nwp], {1:nw, []})';
            E_perf = screen_system([nw, nwp], {[], 1:nwp})';

            %output indexers
            nzp = diss.spec.nzp;
            nz = ssize(diss.plant.D, 1) - nzp - nwp; %for analysis program, the outputs contain copies of the inputs
            E_out = screen_system([nz, nwp, nzp], {1:nz, [],  1:nzp});
           
            plant_rob = E_out * diss.plant*E_rob;

            if nzp > 0
                plant_perf = E_out * diss.plant*E_perf;            
            else
                plant_perf = [];
            end


        end

        function [con_M_out, objective] = quad_performance_augment(obj, diss, vars, con_M, plant_perf)
            %apply a quadratic performance constraint by Schur-Complement
            %in Analysis
            %
            %Args:            
            %   diss (diss_data):   current dissipation constraint
            %   vars:               variables of the problem
            %   con_M (lmim):       current PSD constraint (LMI)
            %   plant_perf (sdpss):  plant with performance outputs
            %
            %Returns:
            %   con_M_out (lmim):   Schur-complemented LMI constraint with
            %                       robustness and performance
            %   objective:          the objective in optimization (scalar)


            if diss.spec.nwp + diss.spec.nzp > 0                
                %an extra quadratic performance condition is present

                %get the supply rate
                vars_spec = vars.spec{diss.spec.id};
                [quad_perf, objective] = diss.spec.supply_quad(vars_spec);                
                nt = ssize(quad_perf.U, 1);

                %apply the schur-complement-style constraint
                perfb = obj.perf_block(plant_perf, quad_perf);

                con_M_out = blkdiag(con_M, zeros(nt));
                con_M_out = con_M_out- perfb;
            else
                %no extra quadratic performance condition 
                objective = 0;
                con_M_out = con_M;
                
            end

        end

    end

    methods (Abstract)
        %variable creation routines                              
    end
end

