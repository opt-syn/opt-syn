classdef op_sml < op_sml_interface
    %OP_SML An operator which is the subdifferential of a function in SmL:
    %
    %F = partial f, where f(x) - m norm(x, 2)^2 and L norm(x, 2)^2 - f(x)
    %are both proper, convex, and closed with -Inf < m <= L < inf
    %
    %
    % noncausal multipliers    
    methods
        function obj = op_sml(m, L, id)
            %OP_SML Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 3
                id = 0;
            end
            obj@op_sml_interface(m ,L, id)      

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
            Cf1 = lmim(['Cf1_', obj.sid], reps, order(1)*reps, 'full');
            Cf2 = lmim(['Cf2_', obj.sid], reps, order(2)*reps, 'full');
            Df1 = lmim(['Df1_', obj.sid], reps, reps, 'full');           
            Df2 = lmim(['Df2_', obj.sid], reps, reps, 'full');           

            
            % orep = order(2)*reps;
            E = lmim(['E_', obj.sid], order(2)*reps, order(1)*reps, 'full');           

            vars = struct('Cf1', Cf1, 'Cf2', Cf2, 'Df1', Df1, 'Df2', Df2, 'E', E);

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

        function cons = filter_constraints(obj, cons, order, vars, iqc)
            %FILTER_CONSTRAINTS constraints on the filter coefficients            

            %Zames-Falb DHD constraints with terminal cost

            h = sum(order) + 1;

            %get the lifted system
            reps = dim(iqc.Psi1.D, 2);
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


            [cons] = dhd_impose_half(P, cons, obj.LMILAB);
            [cons] = dhd_impose_half(P', cons, obj.LMILAB);
        end

        function M = build_M(obj, vars, order, reps);
            %BUILD_M create the running cost M
            M = kron([0, 0, 0, 1; ...
                      0, 0, 1, 0; ...
                      0, 1, 0, 0; ...
                      1, 0, 0, 0], eye(reps));
        end

        function X_out = build_X(obj, vars, order, reps)
            %BUILD_X create the terminal cost X

            %X = [0, E; E', 0];

            [nE, mE] = dim(vars.E);

            znE = zeros(nE);
            zmE = zeros(mE);
            X_out = [znE, vars.E; vars.E', zmE];
        end

        

    end
end

