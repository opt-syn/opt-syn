classdef op_sml < op_sml_interface
    %OP_SML An operator which is the subdifferential of a function in SmL:
    %
    %F = partial f, where f(x) - m norm(x, 2)^2 and L norm(x, 2)^2 - f(x)
    %are both proper, convex, and closed with -Inf < m <= L < inf
    %
    %
    % noncausal multipliers    



    methods
        function obj = op_sml(m, L, c)
            %OP_SML Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 3
                c = 1;
            end
            obj@op_sml_interface(m ,L, c)      

        end
        
        function [vars] = create_vars(obj, order, reps)
            %CREATE_VARS form the variables in an IQC
            %
            %Input: 
            %   order:  order of the IQC [causal, noncausal]
            %   rep:    number of repetitions of the operator
 
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
            %a normalization term for the multipliers
            % cs = trace(vars.Df1) + trace(vars.Df2);
            if obj.same || obj.ERGODIC
                cs = 1;
            else
                cs = trace(vars.Df1);
            end
        end


        function M_erg = ergodic_supply(obj, reps)
            %ERGODIC_SUPPLY supply rate for function value decrease
            %ergodic convergence

            if nargin < 2
                reps = 1;
            end

            M_erg = ergodic_supply@op_sml_interface(reps);

            %TODO: finish and verify this?
        end

        function [Psi1, Psi2] = build_psi(obj, vars, order, reps)
            %BUILD_PSI construct the filter for the SML function
            %
            %use Positive-Real multipliers to do this

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

        function sm = same(obj)
            %SAME: is there any uncertainty in this oracle?
            sm = (obj.m == obj.L);
        end

        function sm = get_same(obj, reps)
            %GET_SAME: is there any uncertainty in this oracle?
            sm = kron(obj.m, eye(reps));
        end


        function [Psi1, Psi2] = build_psi_reduced(obj, vars, order, reps)
            %BUILD_PSI construct the filter for the SML function
            %
            %use Positive-Real multipliers to do this

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
            %get the lifted system for the DHD expression in Zames-Falb
            %guarantees
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

            M0 = [0, 0, 0, 1; ...
                  0, 0, 1, 0; ...
                  0, 1, 0, 0; ...
                  1, 0, 0, 0];
            if obj.SUBLINEAR && ~obj.same
                % Msub = 0;

                sig = 1/(obj.L - obj.m);
                Msub0 = obj.m * [1, sig; sig, sig^2] + sig*[0, 0; 0, 1];
                
                I0rep = diag([1, zeros(reps-1, 1)]);
                Msub = kron(Msub0, I0rep);

                Msub = circshift(kron(Msub, eye(2)), reps, 2);
            else
                Msub = zeros(4*reps);
            end

            M = kron(M0, eye(reps)) + Msub;
        end

        function X_out = build_X(obj, vars, order, reps)
            %BUILD_X create the terminal cost X

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

        

    end
end

