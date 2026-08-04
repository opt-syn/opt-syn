classdef op_quad_causal < op_sml_causal
    %OP_QUAD_CAUSAL a gradient of a quadratic function 1/2 x' Q x, with eigenvalues of Q between
    %  m and L.
    
 
    
    methods
        function obj = op_quad_causal(m, L, c)
            %OP_QUAD_CAUSAL constructor
            

            if nargin < 3
                c = 1;
            end
            obj@op_sml_causal(m , L, c)            

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
                order = 0;
            end
            
            if nargin < 3
                reps = 1;
            end


            %filter coefficients
            Cf = lmim(['Cf_', obj.sid], reps, order*reps, 'full');
            Df = lmim(['Df_', obj.sid], reps, reps, 'full');

            %certificate of positive-realness of filter
            Pf = lmim(['Pf_', obj.sid], order*reps, order*reps, 'sym');


            vars = struct('Cf', Cf, 'Df', Df, 'Pf', Pf);

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
            if isempty(vars)
                cs = 0;
            else                
                cs = trace(vars.Df);
            end
        end

        function cons = filter_constraints(obj, cons, order, vars, iqc)
            %FILTER_CONSTRAINTS constraints on the filter coefficients            
            %positive-real constraints with terminal cost
            %
            %Args:
            %   cons:   accumulated constraints
            %   vars:   variables of the problem             
            %   rho_sched:  which times should be discounted
            %   iqc_out:    the IQC under consideration            
            %Returns:
            %   cons:   accumulated constraints

            [n, m] = size(iqc.Psi1.B);
            Ablock = [eye(n), zeros(n, m);
            iqc.Psi1.A, iqc.Psi1.B];
            Mblock = blkdiag(vars.Pf, -vars.Pf);

            [p, ~] = dim(iqc.Psi1.D);
            sp = [zeros(n, n+p);
            (-0.5)*iqc.Psi1.C, (-0.5)*(iqc.Psi1.D)];

            supply_block = sp + sp';
            % supply_block = sp*2;
            sys_block = Ablock'*Mblock*Ablock;
            
            psd_block = -supply_block + sys_block;

            lmi_curr = [psd_block];

            cons = append_lmi(cons, lmi_curr, obj.LMILAB);

        end

        function [Psi1, Psi2] = build_psi(obj, vars, order, reps)
            %BUILD_PSI construct the filter for the SML function
            %
            %Args:
            %   vars:   variables of the problem    
            %   order:  order of the IQC [number of lags]
            %   reps:    number of repetitions of the operator (from the bind)
            %
            %Returns:
            %   psi1: filter on output (causal)
            %   psi2: filter on input (noncausal components)

            [Af0, Bf0] = block_fir(order);
            Af = kron(eye(reps), Af0 );
            Bf = kron(eye(reps), Bf0);                      
            Cf = vars.Cf;
            Df = vars.Df;
            
            Psi1 = sdpss(Af, Bf, Cf, Df);
            Psi2 = ss(eye(reps));
        end
    end
end

