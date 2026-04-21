classdef op_quad < op_sml
    %OP_QUAD quadratic function 1/2 x' Q x, with eigenvalues of Q between
    %  m and L.
    
 
    
    methods
        function obj = op_quad(m, L, id)
            %OP_QUAD Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 3
                id = 0;
            end
            obj@op_sml(m , L, id)            

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
        
        function cons = filter_constraints(obj, cons, order, vars, iqc)
            %FILTER_CONSTRAINTS constraints on the filter coefficients            

            %Positive-Real constraints with terminal cost
            
            %TODO: validate that these constraints are correct
            reps = dim(iqc.Psi1.D, 2);
            [Psi1, Psi2] = build_psi_reduced(obj, vars, order, reps);
            Psi = blkdiag(Psi1, Psi2);

            [n, m] = size(Psi.B);
            Ablock = [eye(n), zeros(n, m);
            Psi.A, Psi.B];
            Mblock = blkdiag(vars.Pf, -vars.Pf);

            [p, ~] = dim(Psi.D);
            sp = [zeros(n), (-0.5)*Psi.C';
            (-0.5)*Psi.C, (-0.5)*(Psi.D + Psi.D')];

            % supply_block = sp + sp';
            supply_block = sp;
            % supply_block = sp*2;
            sys_block = Ablock'*Mblock*Ablock;
            
            psd_block = -supply_block + sys_block;

            lmi_pass = [psd_block];
            lmi_terminal = vars.Pf - iqc.X;

            cons = append_lmi(cons, lmi_pass, obj.LMILAB);
            cons = append_lmi(cons, lmi_terminal, obj.LMILAB);

            

        end

    end
end

