classdef opt_analysis < opt_manager_interface
    %OPT_ANALYSIS  analysis of optimization algorithms
    %
    % iterative procedure to find a point beta satisfying
    % the fixed-point equation 
    %               0 \in sum_i F_i(\beta).
    %
    properties
        task = 'analysis';
    end
       
    methods
        function obj = opt_analysis(sys)
            %OPT_ANALYSIS Construct an instance of this class            
           
            obj@opt_manager_interface(sys);

            obj.lmi = obj.select_lmi(sys);

        end
        
        %% define IQCs for the operators
        function [obj, vars, cons] = oracle_order(obj,order, ind)
            %ORACLE_ORDER: set the orders of the IQCs
            %Example: order 3 for monotone operators (op_gen) or for SmL
            %causal (op_sml_causal)
            %
            %         order [2, 1] for SmL noncausal (op_sml)

            nop = length(obj.sys.op);
            if ~iscell(order)
                order0 =order;
                order = cell(nop, 1);
                for i = 1:nop
                    order{i}  = order0;
                end
            end

            
            if nargin < 3
                ind = 1:nop;
            end
            vars = obj.vars;
            cons = obj.cons;
            
            for i = 1:nop
                if ismember(i, ind)
                    rep_curr = nnz(obj.sys.bind==i);
                    [iqc_curr, vars_curr,cons_curr] = obj.sys.op{i}.create_iqc(cons, order{i}, rep_curr);
    
                    obj.iqc_op{i} = iqc_curr;
                    vars.op{i} = vars_curr;
                    cons = cons_curr;
                end
            end

            %normalize the coefficients for the filters
            cons = obj.coeff_normalize(vars, cons);

            %TODO: semi-global interface? not very nice.
            obj.cons = cons;
            obj.vars = vars;
        end


        function cons = coeff_normalize(obj, vars, cons)
            %COEFF_NORMALIZE add constraint to normalize the psi multipliers
            
            nop = length(obj.sys.op);
            cs = 0;
            for i = 1:nop
                cs_curr = obj.sys.op{i}.csum_psi(vars.op{i});
                cs = cs + cs_curr;
            end

            cons = append_lmi(cons, cs - nop*0.9, obj.LMILAB);
            cons = append_lmi(cons, -cs + nop*1.1, obj.LMILAB);
        end



        

        %% extract the solution                   
        function  sol = process_recovery(obj, sol, lmi_out);
            %PROCESS_RECOVERY recover the IQCs from the solution
            %
            
            iqc_rec = cell(size(obj.iqc_op));
            for i = 1:length(obj.iqc_op)
                if isnumeric(obj.iqc_op{i})
                    %the Same oracle (m=L, known linear transformation)
                    iqc_rec{i} = obj.iqc_op{i};
                else
                    iqc_rec{i} = obj.iqc_op{i}.recover(lmi_out);
                end

            end

            sol.iqc = iqc_rec;
        end
        

    
    end
end

