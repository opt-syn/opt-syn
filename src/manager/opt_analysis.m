classdef opt_analysis < opt_manager_interface
    %OPT_ANALYSIS  analysis of optimization algorithms
    %
    % iterative procedure to find a point :math:`\beta^*` satisfying
    % the fixed-point equation :math:`0 \in \sum_{i=1}^s F_i(\beta^*)`.

    properties
        order = []; %order of the IQCs
        schedule = []; %when to discount (optionally used for switching)
    end

    methods
        function obj = opt_analysis(sys, config)
            %OPT_ANALYSIS Constructor for analysis
            %
            % Args:
            %   sys: algorithmic system
            %   config: configuration options
                 
            if nargin < 2
                config = opt_config;
            end

            obj@opt_manager_interface(sys, config);

            obj.task = 'analysis';
            obj.lmi = obj.select_lmi(sys);
    
        end
        
        %% define IQCs for the operators
        function obj = process_argument(obj, order)
            %PROCESS_ARGUMENT assign orders to the operators/IQCs
            %
            % Args:
            %   order (cell): orders of the iqcs to search over            
                        
            
            if ~iscell(order)
                %declare the orders
                order0 =order;
                nop = length(obj.sys.op);
                order = cell(nop, 1);
                for i = 1:nop
                    order{i}  = order0;
                end

                
            end


            %check if any are noncausal
            ANY_NONCAUSAL = false;
            for i = 1:length(obj.sys.op)
                if (isa(obj.sys.op{i}, 'op_gen') && (order{i}>0)) || (length(order(i))==2 && (order(2) > 0))
                    ANY_NONCAUSAL = true;
                end
            end
            % % override to enforce the same rho in all channels
            % if ANY_NONCAUSAL
            %     obj.config.gen.same_rho = true;
            %     obj.lmi.config.gen.same_rho = true;
            % end

            obj = obj.oracle_order(order);
        end



        function [obj, vars, cons] = oracle_order(obj,order, ind)
            %ORACLE_ORDER create IQCs at the specified orders
            %
            % Args:
            %   order (cell): orders of the iqcs to search over     
            %   ind (int): which operator to create 
            %
            % Returns
            %   obj: object
            %   vars:   variables of the problem        
            %   cons:       accumulated constraints
            
            
            %Example: order 3 for monotone operators (op_gen) or for SmL
            %causal (op_sml_causal)            
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
            

            %get the discounting schedule (exponents of rho)
            omax = max(cellfun(@(c) sum(c)+1, order));

            obj.schedule = obj.sys.discount_schedule(omax);

            %TODO: semi-global interface? not very nice.
            obj.cons = cons;
            obj.vars = vars;            
            obj.order = order;
            
        end


        function cons = coeff_normalize(obj, vars, cons)
            %COEFF_NORMALIZE add constraint to normalize the psi multipliers
            %
            %Args:
            %   vars:   variables of the problem        
            %   cons:   accumulated constraints
            %
            %Returns:            
            %   cons:    accumulated constraints
            nop = length(obj.sys.op);
            cs = 0;
            for i = 1:nop
                cs_curr = obj.sys.op{i}.csum_psi(vars.op{i});
                cs = cs + cs_curr;
            end


            %LMILAB doesn't like equality constraints on coefficients
            %so add a margin
            marg = obj.config.ana.normalize_margin;

            cons = append_lmi(cons, cs - nop*(1-marg), obj.config.LMILAB);
            cons = append_lmi(cons, -cs + nop*(1+marg), obj.config.LMILAB);
        end



        function [vars, cons, objective, alg_psi, diss] = build_program(obj, specs)
            %form the analysis program
            %Args:
            %   specs (cell): performance specifications
            %Returns:
            %   vars:   variables of the problem        
            %   cons:       accumulated constraints
            %   objective:  single value to be minimized in inner loop (not
            %               the outer loop of bisection)            
            %   alg_psi:  generalized plant                 
            %   diss (diss_data):   dissipation constraints


            if nargin < 2
                specs = obj.specs;
            end

            %BUILD_PROGRAM set up the algorithm analysis problem
            [vars, cons, objective, alg_psi, diss] = build_program@opt_manager_interface(obj, specs); 


            %load in the filter constraints      

            %this requires a weighting by the exponential discounts            
            
            % rho_common = obj.get_common_rho(specs);
            rho_common = 1;
            rho_pow = rho_common.^(obj.schedule);
        
            
            %see if this can be parameterized later
            for i = 1:length(obj.sys.op)
                cons = obj.sys.op{i}.filter_constraints(cons, obj.order{i}, vars.op{i}, rho_pow, obj.iqc_op{i});
            end


            if isscalar(specs) && isa(specs{1}, 'spec_stability')
                cons = obj.coeff_normalize(vars, cons);
            end

        end
           
        function [diss] = index_specs(obj, alg_psi, iqc_data, specs)
            %INDEX_SPECS  index into the performance specifications and
            %form a dissipation relation
            %                       
            %Args:
            %   alg_psi:  generalized plant
            %   iqc_data:  container for the iqcs
            %   specs (cell): performance specifications
            %Returns:                 
            %   diss (diss_data):   dissipation constraints

            
            %TODO: maybe this should go inside the (system), not (manager)?
            
            if nargin < 4
                specs = obj.specs;
            end

            if nargin < 5
                target_ind = 0;
            end

            iqc_op = iqc_data.iqc_op;

            diss = cell(length(specs), 1);
            %determine the indices for each performance specification
            for i = 1:length(specs)
                
                      
                sp = specs{i};
 

                %output indexer                
                if iscell(alg_psi)
                    %TODO: change to genplant_poly type?
                    nz = alg_psi{1}.nz;
                    nw = alg_psi{1}.nw;                    
                    nwp = alg_psi{1}.nwp;
                    nzp = alg_psi{1}.nzp;
                    nu = 0;
                    ny = 0;
                else
                    nz = alg_psi.nz;
                    nw = alg_psi.nw;
                    nwp = alg_psi.nwp;
                    nzp = alg_psi.nzp;
                    nu = 0;
                    ny = 0;
                end
                nzpa = length(sp.izp);
                nwpa = length(sp.iwp);
                i_output = 1:(nz+nzpa+ny);
                j_output = [(1:nz), (nz + sp.izp), nz+nzp + (1:ny)];
                v_output = ones(length(i_output), 1);
                E_output = full(sparse(i_output, j_output, v_output, nz+nzpa+ny, nz+nzp+ny));
                
                %input indexer
                i_input = 1:(nw+nwpa+nu);
                j_input = [(1:nw), (nw + sp.iwp), nw+nwp + (1:nu)];
                v_input = ones(length(i_input), 1);
                E_input = full(sparse(i_input, j_input, v_input, nw+nwpa+nu, nw+nwp+nu))';
                
                %nonminimal representation of the multiplier-extended plant
                if iscell(alg_psi)

                    alg_screen = cell(size(alg_psi));
                    for j = 1:length(alg_screen)
                        alg_screen{j} = E_output* alg_psi.P{j} * E_input;
                    end
                else
                    alg_screen = E_output * alg_psi.P * E_input;
                end


                diss{i} = diss_data;                
                diss{i}.iqc_rob = iqc_op;
                diss{i}.spec = sp;                    
                diss{i}.plant = alg_screen;
                diss{i}.rho = sp.rho;
                % %need to permute the entries of Mdiag for the partition
            end

        end



        

        %% extract the solution                   
        function  sol = process_recovery(obj, sol, lmi_out, alg_psi, diss);
            %PROCESS_RECOVERY recover the IQCs from the solution of the
            %analysis program
            %
            %Args:
            %   sol:  solution structure
            %   lmi_out:  output of solver routines
            %   alg_psi: generalized plant
            %   diss (diss_data):   dissipation constraints
            %Returns:                 
            %   sol:  solution structure
            

            
            iqc_rec = cell(size(obj.iqc_op));
            for i = 1:length(obj.iqc_op)
                if isnumeric(obj.iqc_op{i})
                    %the Same oracle (m=L, known linear transformation)
                    iqc_rec{i} = obj.iqc_op{i};
                else
                    iqc_rec{i} = obj.iqc_op{i}.recover(lmi_out);
                end

            end

            %store the system
            sol.sys = obj.sys;
            %confirm regulator equation
            sol.regcl = obj.lmi.reg.check_regulator();

            %store the IQCs and closed-loop certified systems
            sol.cert.iqc_op = iqc_rec;
            alg_psi_rec = alg_psi;

            if ~iscell(alg_psi)
                alg_psi_rec = {alg_psi_rec};
            end
            
            
            for i = 1:length(alg_psi_rec)
                alg_psi_rec{i}.P.C = double(double(alg_psi_rec{i}.P.C, lmi_out));
                alg_psi_rec{i}.P.D = double(double(alg_psi_rec{i}.P.D, lmi_out));
                alg_psi_rec{i} = ss(alg_psi_rec{i}.A, alg_psi_rec{i}.B, alg_psi_rec{i}.C, alg_psi_rec{i}.D, 1);
            end

            if ~iscell(alg_psi)
                sol.cert.alg_psi = alg_psi_rec{1};
            else
                sol.cert.alg_psi = alg_psi_rec;
            end

            
        end
        

    
    end
end

