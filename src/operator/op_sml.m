classdef op_sml < op_sml_interface
    %OP_SML An operator which is the subdifferential of a function in SmL:
    

    %F = partial f, where f(x) - m norm(x, 2)^2 and L norm(x, 2)^2 - f(x)
    %are both proper, convex, and closed with -Inf < m <= L < inf
    %
    %
    % noncausal multipliers    



    methods
        function obj = op_sml(m, L, c)
            %OP_SML Constructor
            if nargin < 3
                c = 1;
            end
            obj@op_sml_interface(m ,L, c)      

        end
        
        function [vars] = create_vars(obj, order, reps)
            %CREATE_VARS form the variables in an IQC
            %
            %
            %Args: 
            %   order:  order of the IQC [number of lags]
            %   reps:    number of repetitions of the operator (from the bind)
            %
            %Returns:
            %   vars:   variables of the problem            


            if nargin < 2
                order = [0, 0];
            end

            if length(order)==1
                order = [order, 0];
            end
            
            if nargin < 3
                reps = 1;
            end


            %filter coefficients
            if order(1)
                Cf1 = lmim(['Cf1_', obj.sid], reps, order(1)*reps, 'full');
            else
                Cf1 = [];
            end

            if order(2)
                Cf2 = lmim(['Cf2_', obj.sid], reps, order(2)*reps, 'full');
            else
                Cf2 = [];
            end

            Df1 = lmim(['Df1_', obj.sid], reps, reps, 'full');    


            %TODO: allow for general Df2?
            % Df2 = lmim(['Df2_', obj.sid], reps, reps, 'full');           
            Df2 = zeros(reps);

            
            % orep = order(2)*reps;
            if (order(1) > 0) && (order(2) > 0)
                E = lmim(['E_', obj.sid], order(2)*reps, order(1)*reps, 'full');           
            else
                %get the right empty dimensions
                E = zeros(order(2)*reps, order(1)*reps);
            end

            vars = struct('Cf1', Cf1, 'Cf2', Cf2, 'Df1', Df1, 'Df2', Df2, 'E', E);

        end  
        

        function cs = csum_psi(obj, vars)
            %a normalization term for the coefficients, reducing degrees         
            %of freedom in the Analysis problem
            %
            %Args:
            %   vars:   variables of the problem 
            %Returns:
            %   cs: the sum of nonnegative variables
            %

            if obj.same || obj.ERGODIC
                cs = 1;
            else
                cs = trace(vars.Df1);
            end
        end

        function [Psi1, Psi2] = build_psi(obj, vars, order, reps)
            %BUILD_PSI construct the zames-falb filter for the SML function
            %
            %
            %Args:
            %   vars:   variables of the problem    
            %   order:  order of the IQC [number of lags]
            %   reps:    number of repetitions of the operator (from the bind)
            %
            %Returns:
            %   psi1: filter on output (causal)
            %   psi2: filter on input (noncausal components)


            %primal filter
            [Af10, Bf10] = block_fir(order(1));
            Af1 = kron(eye(reps), Af10 );
            Bf1 = kron(eye(reps), Bf10);                      
            Cf1 = [vars.Cf1;zeros(reps, order(1)*reps)];
            Df1 = [vars.Df1; eye(reps)];
            
            %dual filter
            [Af20, Bf20] = block_fir(order(2));
            Af2 = kron(eye(reps), Af20 );
            Bf2 = kron(eye(reps), Bf20);                      
            Cf2 = [vars.Cf2;zeros(reps, order(2)*reps)];
            Df2 = [vars.Df2; eye(reps)];


            Psi1 = sdpss(Af1, Bf1, Cf1, Df1);
            Psi2 = sdpss(Af2, Bf2, Cf2, Df2);    

        end


        function [Psi1, Psi2] = build_psi_reduced(obj, vars, order, reps)
            %BUILD_PSI_REDUCED construct the filter for the SML function
            %but without the identity term (second channel)
            %
            %Args:
            %   vars:   variables of the problem    
            %   order:  order of the IQC [number of lags]
            %   reps:    number of repetitions of the operator (from the bind)
            %
            %Returns:
            %   psi1: filter on output (causal)
            %   psi2: filter on input (noncausal components)

            %primal filter
            [Af10, Bf10] = block_fir(order(1));
            Af1 = kron(eye(reps), Af10 );
            Bf1 = kron(eye(reps), Bf10);                      
            Cf1 = [vars.Cf1 ];
            Df1 = [vars.Df1 ];
            
            %dual filter
            [Af20, Bf20] = block_fir(order(2));
            Af2 = kron(eye(reps), Af20 );
            Bf2 = kron(eye(reps), Bf20);                      
            Cf2 = [vars.Cf2 ];
            Df2 = [vars.Df2 ];


            Psi1 = sdpss(Af1, Bf1, Cf1, Df1);
            Psi2 = sdpss(Af2, Bf2, Cf2, Df2);    

        end

        function P = dhd_lift(obj, order, vars, iqc)
            %DHD_LIFT get the lifted system for the doubly hyperdominant 
            % expression in Zames-Falb guarantees
            %
            %Args:            
            %   order:  order of the IQC [number of lags]
            %   vars:   variables of the problem    
            %   iqc:  the iqc under consideration            
            %
            %Returns:
            %   P:  matrix that should be DHD 

            if ~isnumeric(iqc.Psi1.D)
                reps = dim(iqc.Psi1.D, 2);
            else
                reps = size(iqc.Psi1.D, 2);
            end
            h = sum(order) + 1;

            [Psi1_red, Psi2_red] = build_psi_reduced(obj, vars, order, reps);
                        
            Psi1_lift = sdpsslift(Psi1_red, h);
            Psi2_lift = sdpsslift(Psi2_red, h);

            Dh1 = Psi1_lift.D;
            Dh2 = Psi2_lift.D;

            Bh1 = Psi1_lift.B;
            Bh2 = Psi2_lift.B;

            %form the DHD constraint
            E = vars.E;

            BE = Bh2'*E*Bh1;

            P = Dh1 + Dh2' + BE;      
        end

        function cons = filter_constraints(obj, cons, order, vars, rho_sched, iqc)
            %FILTER_CONSTRAINTS constraints on the filter coefficients            
            %Zames-Falb DHD constraints with terminal cost
            %
            %Args:
            %   cons:   accumulated constraints
            %   vars:   variables of the problem             
            %   rho_sched:  which times should be discounted
            %   iqc_out:    the IQC under consideration            
            %Returns:
            %   cons:   accumulated constraints


            if ~isscalar(iqc)
                if isscalar(order)
                    order = [order, 0];
                end
    
                reps = size(iqc.Psi1.B, 2)/obj.c;
                
                    
                P = obj.dhd_lift(order, vars, iqc);
    
                %impose the exponential discounting
                nsched = size(rho_sched, 2);
                h = sum(order)+1;
                for i = 1:nsched
                    % rho_1 = kron(diag(rho_sched(1:(order(1)+1), i)), eye(reps));
                    % rho_2 = kron(diag(rho_sched(1:(order(2)+1), i)), eye(reps));
                    rho_1 = kron(diag(rho_sched(1:h, i)), eye(reps));
                    rho_2 = rho_1;
    
                    P_rho = rho_1 * P * rho_2;
    
                    [cons] = dhd_impose(P_rho, cons, obj.LMILAB);
                end
                
                nu = dim(vars.Df1, 1);
                for i = 1:nu
                    %TODO: normalization of the multipliers
                    e = zeros(nu, 1);
                    e(i) = 1;
                    cons = append_lmi(cons, e'*vars.Df1*e, obj.LMILAB);
                    % cons = append_lmi(cons, e'*vars.Df2*e, obj.LMILAB);
                end
            end
        end

        function M = build_M(obj, vars, order, reps);
            %BUILD_M create the running cost M
            %Args:
            %   vars:   variables of the problem 
            %   order:  order of the IQC [number of lags]
            %   reps:    number of repetitions of the operator (from the bind)
            %Returns:
            %   M_out: the running cost
        
            M0 = [0, 0, 0, 1; ...
                  0, 0, 1, 0; ...
                  0, 1, 0, 0; ...
                  1, 0, 0, 0];
            if obj.ERGODIC && ~obj.same
                % Msub = 0;

                Msub0 = obj.ergodic_supply(reps);

                I0rep = diag([0, 1]);
                Msub = kron(Msub0, I0rep);

                % Msub = kron(Msub, eye(2));
            else
                Msub = zeros(4*reps);
            end

            M = kron(M0, eye(reps)) + Msub;
        end

        function X_out = build_X(obj, vars, order, reps)
            %BUILD_X create the terminal cost X
            %
            %Args:
            %   vars:   variables of the problem 
            %   order:  order of the IQC [number of lags]
            %   reps:    number of repetitions of the operator (from the bind)
            %Returns:
            %   X_out: the terminal cost


            %X = [0, E; E', 0];

            
            if isempty(vars.E)
                [nE, mE] = size(vars.E);
                X_out = zeros(max(nE, mE));
            else
                [nE, mE] = dim(vars.E);
    
                znE = zeros(nE);
                zmE = zeros(mE);
                X_out = [znE, vars.E; vars.E', zmE];
            end
        end

        function [iqc] = create_iqc_identity(obj, reps)
            %CREATE_IQC_IDENTITY form a valid IQC satisfied by the sml
            %operator. This is used as a warm start in synthesis.
            %
            %Args:             
            %   reps:    number of repetitions of the operator (from the bind)
            %
            %Returns:
            %   iqc (iqc_loop_split): a valid IQC with no dynamics    

            iqc = create_iqc_identity@op_sml_interface(obj, reps);

            if ~isnumeric(iqc) && obj.ERGODIC

                Msub = obj.ergodic_supply(reps);

                iqc.M = iqc.M + Msub;

            end
        

        end

        

    end
end

