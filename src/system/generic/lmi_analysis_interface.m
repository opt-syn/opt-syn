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
        function obj = lmi_analysis_interface(sys)
            %LMI_ANALYSIS_INTERFACE Construct an instance of this class
            %   Detailed explanation goes here
            obj@lmi_dispatch_interface(sys);
        end


        %% variable creation

        function [vars, cons] = create_vars(obj, vars, cons, alg_psi, specs)
            %CREATE_VARS create the variables for the problem

            [vars.diss, cons] = obj.create_vars_storage(cons, alg_psi);
            [vars.spec, cons] = obj.create_vars_spec(cons, specs);
        end




        function [G, cons] = define_storage_G(obj, cons, alg_psi,  name)
            %DEFINE_STORAGE_G storage function for a specific subsystem
            n = ssize(alg_psi{1}.A, 1);
            G = lmim(['G', name], n, n, 'sym');

            

            if obj.tol.G_max < Inf    
                %issue in the bounding?
                cons = append_lmi(cons, obj.tol.G_max*eye(n)  - G, obj.LMILAB);
                cons = append_lmi(cons, obj.tol.G_max*eye(n)  + G, obj.LMILAB);                
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
                [vars_spec{i}, cons] = specs{i}.create_vars(cons);
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
            cons = append_lmi(cons, con_X - eye(sx)*obj.tol.X, obj.LMILAB);

        end

        %% helper functions to construct LMIs

        function sb = sys_block(obj, plant, Pnew, Pold, rho)
            % SYS_BLOCK system block used in analysis programs
            %
            %sb =  [0, I]^T [Pold*rho^2, 0] [0, I]
            %      [A, B]   [0,      -Pnew] [A, B]

            
            A = plant.A;
            B = plant.B;
            
            [n, m] = size(B);  

            Ablock = [eye(n), zeros(n, m);
                A, B];

            Pblock = blkdiag(Pold*(rho^2), -Pnew);

            sb = Ablock' * Pblock * Ablock;            

        end

        function sb = supply_block(obj, plant, M)
            % SUPPLY_BLOCK supply block used in analysis programs
            %
            %sb =  [C, D]^T [-M] [C D]
            %               


            Cblock = [plant.C, plant.D];   
            
            sb = Cblock' * (-M) * Cblock;

        end


        function [plant_no_p, CDp] = separate_performance_output(obj, diss)

            %SEPARATE_PERFORMANCE_OUTPUT
            %extract the performance output from the plant
            %
            %plant_no_p:    the plant with the performance output removed
            %CDp:           the entries of [Cp, Dp] matrix for the
            %               performance output
            %
            %
            %This routine is used in the computation of l2 norms via Schur
            %complements            
            
            %separate the performance outputs   
            nzp = length(diss.spec.izp);
            ind_sep = (diss.iqc_rob.np) + (1:nzp);
            nz = ssize(diss.plant.D, 1);

            ind_diff_sep = setdiff(1:nz, ind_sep);
            Iz = eye(nz);
            Izp = Iz(ind_diff_sep, :);

            %the plant without the schur-complemented-out performance input
            plant_no_p = Izp*diss.plant;


            Ezp = sparse(1:length(ind_sep), ind_sep, ones(length(ind_sep)), ...
                length(ind_sep), ssize(diss.plant.C, 1));


            CDp = Ezp * [diss.plant.C, diss.plant.D];

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

