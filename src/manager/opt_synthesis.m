classdef opt_synthesis < opt_manager_interface
    %OPT_SYNTHESIS synthesis of optimization algorithms
    %
    % iterative procedure to find a point :math:`\beta^*` satisfying
    % the fixed-point equation :math:`0 \in \sum_{i=1}^s F_i(\beta^*)`,    
    % in which the oracles :math:`F_i` are interfaced over a dynamical network
    
    methods
        function obj = opt_synthesis(sys, config, iqc_op)
            %OPT_SYNTHESIS Constructor for synthesis
            % Args:
            %   sys: algorithmic system
            %   config: configuration options
            %   iqc_op (cell): of IQCs for the operators
                 

            if nargin < 2
                config = opt_config;
            end

            obj@opt_manager_interface(sys, config);
                        
            if nargin <= 2 || isempty(iqc_op)                
                %create identity filters, but keep the loop transformations
                obj.iqc_op = obj.make_blank_iqc;
            else
                obj.iqc_op = iqc_op;            
            end

            obj.task = 'synthesis';            
            obj.lmi = obj.select_lmi(sys);

        end

        function iqc_op = make_blank_iqc(obj)
            %if no IQCs are provided, make identity IQCs
            %
            % Returns:
            %   iqc_op (cell): IQCs for the operators
                nop = length(obj.sys.op);
                iqc_op = cell(nop, 1);
                bind =obj.sys.bind;
                for i = 1:nop
                    nrep = sum(i==bind);
                    op_blank = obj.sys.op{i}.create_iqc_identity(nrep);
                    
                    if ~isnumeric(op_blank)
                        iqc_op{i} = op_blank.factor();
                    else
                        iqc_op{i} = op_blank;
                    end
                end
        end
        
        function obj = process_argument(obj,iqc_op)
            %PROCESS_ARGUMENT assign orders to the operators/IQCs
            % Args:
            %   iqc_rob: IQCs representing the robust uncertainties

            if nargin > 1 && ~isempty(iqc_op)
                obj.iqc_op = iqc_op;
            else
                obj.iqc_op = obj.make_blank_iqc();
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

            % [rho, sperf] = obj.perf_specs(specs);
            sperf = specs;
            rho = sperf{1}.rho;

            iqc_op = iqc_data.iqc_op;

            diss = cell(length(sperf), 1);
            %determine the indices for each performance specification
            for i = 1:length(sperf)
                
                      
                sp = sperf{i};
                iwp_iqc = (1:(iqc_op.nw))';
                ir_iqc_first = (1:(iqc_op.np))';


                count_iqc_in = (iqc_op.nw);
                count_iqc_out = (iqc_op.np);

                %any performance outputs?
                if isempty(sp.izp)
                    ir_iqc_first_r =[];
                    
                else
                    ir_iqc_first_r = count_iqc_out  + (sp.izp);
                    count_iqc_out = count_iqc_out + length(sp.izp);
                end
                
                %any performance inputs?
                if isempty(sp.iwp)
                    iw_iqc_first_r = [];
                else
                    iw_iqc_first_r = count_iqc_in + (sp.iwp);
                    count_iqc_in = count_iqc_in + length(sp.iwp);                    
                end

                iwp_iqc = [iwp_iqc; iw_iqc_first_r];

                ir_iqc0 = [ir_iqc_first; ir_iqc_first_r];
               
                sp_ind_w = iwp_iqc;
                sp_ind_r = ir_iqc0;
                
                if iscell(alg_psi)
                    %TODO: change to genplant_poly type?
                    nwr = alg_psi{1}.nz;
                    nww = alg_psi{1}.nw;                    
                    nu = alg_psi{1}.nu;
                    ny = alg_psi{1}.ny;
                else
                    nwr = alg_psi.nz;
                    nww = alg_psi.nw;
                    nu = alg_psi.nu;
                    ny = alg_psi.ny;
                end

                nwr = nwr + length(sp.izp);
                nww = nww + length(sp.iwp);

                
                E_r = blkdiag(full(sparse(1:length(sp_ind_r), sp_ind_r, ...
                    ones(1, length(sp_ind_r)), length(sp_ind_r), nwr)), eye(ny));

                E_w = blkdiag(full(sparse(1:length(sp_ind_w), sp_ind_w, ...
                    ones(1, length(sp_ind_w)), length(sp_ind_w), nww)), eye(nu));

                

                %enforce squareness in the performance specs?

                

                %nonminimal representation of the multiplier-extended plant

                %TODO: write fancier index code?
                

                

                if iscell(alg_psi)

                    n2 = alg_psi{1}.dump_dim();
                    n2.nwp = length(sp.iwp);
                    n2.nzp = length(sp.izp);

                    alg_screen = cell(size(alg_psi));
                    for j = 1:length(alg_screen)
                        alg_screen{j} = genplant(E_r * alg_psi{j}.ss * E_w, n2);
                    end                    
                    %TODO: write this part: cells/genplant poly
                else

                    n2 = alg_psi.dump_dim();
                    n2.nwp = length(sp.iwp);
                    n2.nzp = length(sp.izp);

                    alg_screen_P = E_r * alg_psi.ss * E_w;
                    alg_screen = genplant(alg_screen_P, n2);

                    
                    
                end

                % iqc_data_rump = rmfield(iqc_data, 'iqc_op');
                
                diss{i} = diss_data;
                diss{i}.iqc_rob = iqc_op;
                diss{i}.rho = rho;
                diss{i}.spec = sp;
                diss{i}.iqc_data = iqc_data;
                diss{i}.plant = alg_screen;
                diss{i}.plant_reg = obj.lmi.reg.sys_regulated_aug();
                diss{i}.ndiss = length(specs);





                %TODO: this may run into trouble if one entry has an X.
                %performance with dynamic multipliers?
            end

        end

        %% extract the solution                   
        function  sol = process_recovery(obj, sol, lmi_out, alg_psi, diss)
            %PROCESS_RECOVERY recover the IQCs from the solution of the
            %synthesis program
            %
            %Args:
            %   sol:  solution structure
            %   lmi_out:  output of solver routines
            %   alg_psi: generalized plant
            %   diss (diss_data):   dissipation constraints
            %Returns:                 
            %   sol:  solution structure
            sol.cert.iqc_op = obj.iqc_op;
            sol.cert.iqc_op_all = obj.iqc_op_all;
            sol.vars.rho = sol.rho;
            [sol] = obj.lmi.process_recovery(sol, lmi_out, alg_psi, diss);            

            sol.sys = obj.sys;
            sol.sys.K = sol.cert.K;

            %check regulator equation
            reg2 = obj.lmi.reg;
            reg2.sys = sol.sys;
            sol.cert.regcl = reg2.check_regulator();
        end

        %% alternating design
        function [sol_history, vr_history, success] = alternate(obj, Niter, order, iqc_init, specs, b_opts)
            %ALTERNATE alternating synthesis and analysis. 
            % use bisection in analysis and synthesis if rho is minimized.
            %
            %
            %Args:
            %   Niter (int): number of alternation iterations
            %   order (cell):  orders of the operators (for analysis)
            %   iqc_init (cell):  initial IQCs for the operators (for synthesis)            
            %   specs (cell):   performance specifications (for both)            
            %   b_opts:   (bisect_opts) bisection options (bisect_opts)
            %Returns:                 
            %   sol_history (cell):  cell of solutions, first row is Synthesis, second row is analysis.
            %   vr_history (cell):   lower and upper bound of parameter
            %   success (bool):      success of alternation method



            sol_history = cell(2, Niter);
            vr_history = cell(2, Niter);

            if nargin < 4
                iqc_init = [];
            end

            if nargin < 5
                specs = [];
            end

            if nargin < 6
                b_opts = bisect_opts();
            end
            iqc_curr = iqc_init;
            
            sys_curr = obj.sys;


            if ~iscell(specs)
                if isempty(specs)
                    specs = {spec_stability(1)};
                else
                    specs = {specs};
                end
            end

            for j = 1:numel(iqc_curr)
                iqc_curr{j} = iqc_curr{j}.factor();
            end

            success = false;
            for i = 1:Niter
                %start with synthesis
                if b_opts.bisect
                    [sol_syn, vr_syn] = obj.bisect(iqc_curr, specs, b_opts);

                    %back off a bit
                    vr_back = vr_syn(2)+ b_opts.backoff;

                    spec_back = obj.modify_spec(vr_back, specs, b_opts);

                    obj.specs= {};
                    sol_syn_back = obj.solve_single(iqc_curr, spec_back);

                else
                    [sol_syn] = obj.solve_single(iqc_curr, specs);
                    vr_syn = sol_syn.objective * [1, 1];

                    sol_syn_back = sol_syn;
                end

                
                sol_history{1, i} = sol_syn;
                vr_history{1, i} = vr_syn;

                if isempty(sol_syn) || sol_syn.status
                    break
                end
                


                %then do analysis
                sys_curr = sol_syn_back.sys;

                ana = opt_analysis(sys_curr);
                if b_opts.bisect
                    [sol_ana, vr_ana] = ana.bisect(order, specs, b_opts);
                else
                    [sol_ana] = ana.solve_single(order, specs);
                    vr_ana = sol_ana.objective * [1, 1];
                end

                sol_history{2, i} = sol_ana;
                vr_history{2, i} = vr_ana;

                if sol_ana.status
                    break
                elseif i==Niter
                    success = true;
                end
                
                %prepare for next go-around
                %factor the iqcs from analysis for use in synthesis
                iqc_curr = cell(length(sol_ana.cert.iqc_op), 1);
                
                if i < Niter
                    for j = 1:numel(iqc_curr)
                        if isnumeric(sol_ana.cert.iqc_op{j})
                            iqc_curr{j} = sol_ana.cert.iqc_op{j};
                        else
                            iqc_curr{j} = sol_ana.cert.iqc_op{j}.factor();
                        end
                    end        
                end

            end

            


            
        end



        
    end
end

