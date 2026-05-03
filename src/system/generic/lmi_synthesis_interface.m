classdef lmi_synthesis_interface < lmi_dispatch_interface
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
        reg; %internal model of the controller
                
    end

    
    methods
        function obj = lmi_synthesis_interface(sys, config)
            %LMI_SYNTHESIS_INTERFACE Construct an instance of this class
            %   Detailed explanation goes here
            obj@lmi_dispatch_interface(sys, config);


            %form the internal model
            reg_name = ['regulator_', sys.get_type()];

            reg_handle = str2func(reg_name);
            obj.reg = reg_handle(sys);

            %TODO: better options handling down below
            obj.config = config;
        end

        

        %% variable creation

        function [vars, cons] = create_vars(obj, vars, cons, alg_psi, specs)
            %CREATE_VARS create the variables for the problem

            [vars.diss, cons] = obj.create_vars_storage(cons, alg_psi);
            [vars.spec, cons] = obj.create_vars_spec(cons, specs);
            [vars.reg]  = obj.create_vars_regulator();
            [vars.K, cons]    = obj.create_vars_controller(cons, alg_psi);
        end

        function [vars_diss, cons]= create_vars_storage(obj, cons, alg_psi, name)
            %create_vars_storage create variables for the dissipation
            %constraints
            %
            %Input:
            %   cons:       accumulated constraints
            %   alg_psi:    the filtered algorithmic interconnection
            %   name:       a name for the variable

            if nargin < 4
                name = [];
            end

            [GX, GY, cons] = obj.define_storage_G(cons, alg_psi, name);
            vars_diss= struct('GX', GX, 'GY', GY);

        end

        function [vars_reg] = create_vars_regulator(obj)
            %CREATE_VARS_REGULATOR
            %parameterize the solutions to the regulator equations
            %use this as a variable in reduced-order control
            %
            %systems with more outputs than oracles can have freedom in the            
            %regulator equations (such as optimization problems with known 
            % Laplacian matrices)

            

            
            
            
            vars_reg = obj.reg.create_vars();
            

        end



        function [GX, GY, cons] = define_storage_G(obj, cons, alg_psi,  name)
            %DEFINE_STORAGE_G storage function for a specific subsystem

            %without terminal cost:
            %
            %[GX, I;
            %[I, GY] is PD
            %


            n = ssize(alg_psi.A, 1);
            ns = ssize(obj.reg.S, 1);



            nX = n + ns;

            
            if obj.config.syn.reduced_order
                %TODO: not yet implemented
                nY = n;
            else
                nY = n + ns;
            end

             GX = lmim(['GX', name], nX, nX, 'sym');
            
          

            GY = lmim(['GY', name], nY, nY, 'sym');

            %TODO the terminal constraints with the coupling condition?


            %bound the entries of the GX and GY matrices
            if obj.config.tol.G_max < Inf                   
                cons = append_lmi(cons, obj.config.tol.GX_max*eye(nX)  - GX, obj.config.LMILAB);                            
                cons = append_lmi(cons, obj.config.tol.GY_max*eye(nY)  - GY, obj.config.LMILAB);                
            end            
        end

        function G = get_storage(obj, vars_diss, vars_reg)
            %GET_STORAGE get the storage function matrix G

            % 
            %
            GX = vars_diss.GX;
            GY = vars_diss.GY;
            
            nx = ssize(GX, 1);
            % ns = ssize(obj.reg.S, 1);


            if obj.config.syn.reduced_order
                %TODO: not yet implemented
                %some sort of indexing on Pi
                Pi = vars.reg.Pi;
                G = [GY, eye(nx); eye(nx), GX];
            else
                %
                G = [GY, eye(nx); eye(nx), GX];                
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
        
        function [vars_K, cons] = create_vars_controller(obj, cons, alg_psi, name)
            %CREATE_VARS_CONTROLLER create the nonlinearly-transformed
            %controller matrices

            %get the dimensions


            if nargin < 4
                name = [];
            end

            n = ssize(alg_psi.A, 1);
            ns = ssize(obj.reg.S, 1);
            
            if obj.config.syn.reduced_order
                %TODO: not yet implemented
                nc = n;
            else
                nc = n + ns;
            end

            ny = obj.sys.P.ny;
            nu = obj.sys.P.nu;

            %declare the variables
            vars_K = struct;
            %easy: ABC
            vars_K.A = lmim(['Ak', name], nc, nc);
            vars_K.B = lmim(['Bk', name], nc, ny);
            vars_K.C = lmim(['Ck', name], ns + nu, nc);


            vars_K.D = obj.form_Dk(alg_psi);
            %TODO: better interface here: number of inputs
            


            %no elimination just yet

            %bound entries of the controllers
            kq = [vars_K.A, vars_K.B;            
              vars_K.C,  vars_K.D];
            cons= append_lmi(cons, obj.config.tol.K_max*eye(sum(kq.dim)) - [zeros(kq.dim(1)), kq; kq', zeros(kq.dim(2))], obj.LMILAB);

        end

        function [Dk] = form_Dk(obj, alg_psi, name)
            %FORM_Dk: lower triangular structure needed for the controller
            %need a better interface for the mask


            %also, maybe an object structure for the internal model?
            
            if nargin < 3
                name = [];
            end

            
            n = ssize(alg_psi.A, 1);

            s = length(obj.sys.bind);
            c = obj.sys.op{1}.c;
            nu = obj.sys.P.nu;
            ns = size(obj.reg.R, 2);

            %more difficult: Dk
            
            %the unconstrained term for the internal model control
            Dk1_var = lmim(['Dk1', name], ns + n, ns, 'full');
            Dk = Dk1_var;



            %the sparsity-constrained term for internal model control
            % if nargin < 2
                D_mask_0 = obj.config.syn.D_mask;
            % end
           
            if isempty(D_mask_0)
                D_mask_0 = tril(ones(length(obj.sys.bind)));
            end

            

            %use the triangular structure
            
            %TODO: graceful handling of other dynamics
            D_mask = kron(D_mask_0, eye(c));
 
            nd2= nnz(D_mask);
            if nnz(D_mask) > 0
                Dk2_var = lmim(['Dk2', name], 1, nd2, 'full');
    
                %make sure that the Dc2 term of the subcontroller is
                %lower-triangular
                
                counter = 1;
                for i = 1:s*c
                     if any(D_mask(i, :))
                        eind = find(D_mask(i, :));
                        ncc = length(eind);
                        Dvar_mat = sparse(counter + (1:ncc)-1, 1:ncc, ones(ncc, 1), nd2, ncc);

                        Dvar = Dk2_var * Dvar_mat;
                        
                        Din = sparse(1:ncc, eind, ones(ncc, 1), ncc, s);
                        Din_var = Dvar * Din;

                        Dk = [Dk; Din_var];
    
                        
                        % Dk2_curr = sparse(, 1:nc )
                        % vars.Dk(i+(nu-s), j) = Dk2_var * eind;
                        counter = counter+nnz(D_mask(i, :));
                    else
                        Dk = [Dk; zeros(1, ny)];
                    end
                    % end
                    
                % end
                end
            else
                Dk = [Dk; zeros(s*c, s*c)];
            end

        end

        %% terminal constraints        
        function [cons, con_X] = con_terminal(obj, G, cons,  alg_psi, iqc_op)
            %CON_TERMINAL
            %terminal cost constraint (nonnegativity for the storage function G)
            %coupled positivity if the IQC has a terminal cost

            %too many arguments taken here
            X = iqc_op.X;

            %TODO: allow for reduced-order control           

            %TODO: check that this is the right formula, specifically when
            %X is a non-PSD terminal cost

            %matrix dilation results

            nf = ssize(X);
            n = ssize(G, 1);

            Ef = [eye(nf); zeros(n-nf, nf)];


            X_f = Ef * X * Ef';
            con_X = G + X_f;

            sx = ssize(con_X, 1);
            cons = append_lmi(cons, con_X - eye(sx)*obj.config.tol.X, obj.LMILAB);

        end

        %% helper functions to construct LMIs

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

