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
        

       function cons = filter_constraints(obj, cons, order, vars, rho_sched, iqc)
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

            if isscalar(order)
                order = [order, 0];
            end

                
            P = obj.dhd_lift(order, vars, iqc);

            P_sym = P + P';

            cons = append_lmi(cons, P_sym, obj.LMILAB);
        end
        
        % function cons = filter_constraints(obj, cons, order, vars, iqc)
        %     %FILTER_CONSTRAINTS constraints on the filter coefficients            
        % 
        %     %Positive-Real constraints with terminal cost
        % 
        %     %TODO: validate that these constraints are correct
        %     reps = dim(iqc.Psi1.D, 2);
        %     [Psi1, Psi2] = build_psi_reduced(obj, vars, order, reps);
        %     Psi = blkdiag(Psi1, Psi2);
        % 
        %     [n, m] = size(Psi.B);
        %     Ablock = [eye(n), zeros(n, m);
        %     Psi.A, Psi.B];
        %     Mblock = blkdiag(vars.Pf, -vars.Pf);
        % 
        %     [p, ~] = dim(Psi.D);
        %     sp = [zeros(n), (-0.5)*Psi.C';
        %     (-0.5)*Psi.C, (-0.5)*(Psi.D + Psi.D')];
        % 
        %     % supply_block = sp + sp';
        %     supply_block = sp;
        %     % supply_block = sp*2;
        %     sys_block = Ablock'*Mblock*Ablock;
        % 
        %     psd_block = -supply_block + sys_block;
        % 
        %     lmi_pass = [psd_block];
        %     lmi_terminal = vars.Pf - iqc.X;
        % 
        %     cons = append_lmi(cons, lmi_pass, obj.LMILAB);
        %     cons = append_lmi(cons, lmi_terminal, obj.LMILAB);
        % 
        % 
        % 
        % end

    end
end

