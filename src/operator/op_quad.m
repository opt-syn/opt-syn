classdef op_quad < op_sml
    %OP_QUAD a gradient of a quadratic function 1/2 x' Q x, with eigenvalues of Q between
    %  m and L.
        
    methods
        function obj = op_quad(m, L, c)
            %OP_QUAD Constructor
            if nargin < 3
                c = 0;
            end
            obj@op_sml(m , L, c)            

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

            vars = create_vars@op_sml(obj, order, reps);

            Pf = lmim(['Pf_', obj.sid], sum(order)*reps, sum(order)*reps, 'sym');

            vars.Pf = Pf;
        end

        function M = build_M(obj, vars, order, reps);
            %BUILD_M create the running cost M
            %Args:
            %   vars:   variables of the problem 
            %   order:  order of the IQC [number of lags]
            %   reps:    number of repetitions of the operator (from the bind)
            %Returns:
            %   M_out: the running cost

            M0 = kron([0, 1;1, 0], eye(order(1)+1));
            % if obj.ERGODIC && ~obj.same
            %     % Msub = 0;
            % 
            %     Msub0 = obj.ergodic_supply(reps);
            % 
            %     I0rep = diag([0, 1]);
            %     Msub = kron(Msub0, I0rep);
            % 
            %     % Msub = kron(Msub, eye(2));
            % else
            %     Msub = zeros(4*reps);
            % end
            Msub = 0;

            M = kron(eye(reps), M0) + Msub;
        end
        

       % function cons = filter_constraints(obj, cons, order, vars, rho_sched, iqc)
       %      %FILTER_CONSTRAINTS constraints on the filter coefficients                        
       %      %positive-real constraints with terminal cost
       %      %
       %      %Args:
       %      %   cons:   accumulated constraints
       %      %   vars:   variables of the problem             
       %      %   rho_sched:  which times should be discounted
       %      %   iqc_out:    the IQC under consideration            
       %      %Returns:
       %      %   cons:   accumulated constraints
       % 
       %      if isscalar(order)
       %          order = [order, 0];
       %      end
       % 
       % 
       %      P = obj.dhd_lift(order, vars, iqc);
       % 
       %      P_sym = P + P';
       % 
       %      cons = append_lmi(cons, P_sym, obj.LMILAB);
       % end


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


           %this is a convex formulation with extra channels, all
           %coefficients are in 

           %primal filter
           [Af10, Bf10] = block_fir(order(1));
           Af1 = kron(eye(reps), Af10 );
           Bf1 = kron(eye(reps), Bf10);  
           Cf1 = kron(eye(reps), [eye(order(1)); zeros(1, order(1))]);
           Df1 = kron(eye(reps), [zeros(order(1), 1); 1]);
           

           %dual filter
           [Af20, Bf20] = block_fir(order(2));
           Af2 = kron(eye(reps), Af20 );
           Bf2 = kron(eye(reps), Bf20);                      
           Cf2 = [kron(eye(reps), zeros(order(1), order(2))); vars.Cf2];              
           Df2 = [vars.Cf1'; vars.Df1];


           Psi1 = sdpss(Af1, Bf1, Cf1, Df1);
           Psi2 = sdpss(Af2, Bf2, Cf2, Df2);    

       end

       function loop = build_loop(obj, reps)
           %BUILD_LOOP construct the signal transformation matrix
           %
           %Args:
           %   reps:    number of repetitions of the operator (from the bind)
           %
           %Returns:
           %   loop_out: signal transformation matrix for the operator

           if obj.L_top
               loop_base = [-obj.sigma, 1; 1, obj.m];
           else
               loop_base = [-obj.sigma, 1; -1, obj.L];
           end

           loop = kron(loop_base, eye(reps));
       end

       function cons = filter_constraints(obj, cons, order, vars, rho_sched, iqc)
           %FILTER_CONSTRAINTS constraints on the filter coefficients            

           %Positive-Real constraints with terminal cost

           %
           %linear map: no rho weighting needed
           %
           %TODO: validate that these constraints are correct



           %assemble the filter for the sum
           Psi1 = iqc.Psi1;
           Psi2 = iqc.Psi2;
           

           Apsi = blkdiag(Psi1.A, Psi2.A);
           Bpsi = [Psi1.B; Psi2.B];
           Cpsi = blkdiag(Psi1.C, Psi2.C);
           Dpsi = [Psi1.D; Psi2.D];


                    
           %call the KYP lemma for positive realness


           %dynamics
           [n, m] = ssize(Bpsi);
           Ablock = [eye(n), zeros(n, m);
           Apsi, Bpsi];

           Pblock = blkdiag(vars.Pf, -vars.Pf);

           M_sys = Ablock' * Pblock * Ablock;


           %supply
           Cblock = [Cpsi, Dpsi];
           p = ssize(Dpsi, 1);

           sup_block = kron([0, 1; 1, 0], eye(p/2));
           M_supp = Cblock' * sup_block * Cblock;

           %imposition
           lmi_pass = M_sys + M_supp;
           lmi_terminal = vars.Pf - iqc.X;

           
           cons = append_lmi(cons, lmi_pass, obj.LMILAB);
           cons = append_lmi(cons, lmi_terminal, obj.LMILAB);

       end

        
    end
end

