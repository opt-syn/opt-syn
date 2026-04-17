classdef op_quad_causal < op_sml_causal
    %OP_QUAD quadratic function 1/2 x' Q x, with eigenvalues of Q between
    %  m and L.
    
 
    
    methods
        function obj = op_quad_causal(m, L, id)
            %OP_QUAD Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 3
                id = 0;
            end
            obj@op_sml_causal(m , L, id)            

        end
        


        function [vars] = create_vars(obj, order, reps)
            %CREATE_VARS form the variables in an IQC
            %
            %Input: 
            %   order:  order of the IQC [causal, noncausal]
            %   rep:    number of repetitions of the operator
 
            if nargin < 2
                order = 0;
            end
            
            if nargin < 3
                reps = 1;
            end


            %filter coefficients
            Cf = lmim(['Cf', obj.sid], reps, order*reps, 'full');
            Df = lmim(['Df', obj.sid], reps, reps, 'full');

            %certificate of positive-realness of filter
            Pf = lmim(['Pf', obj.sid], order*reps, order*reps, 'sym');


            vars = struct('Cf', Cf, 'Df', Df, 'Pf', Pf);

        end       

        function cons = filter_constraints(obj, cons, vars, iqc)
            %constraints on the filter coefficients

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
            %use Positive-Real multipliers to do this

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

