classdef opt_synthesis < opt_manager_interface
    %OPT_SYNTHESIS synthesis of optimization algorithms
    %
    % iterative procedure to find a point beta satisfying
    % the fixed-point equation 
    %               0 \in sum_i F_i(\beta).
    %
    % in which the oracles F_i are interfaced over a dynamical network
    properties
        % opts = 
        % iqc_op = {};
    end
    
    methods
        function obj = opt_synthesis(sys, iqc_op)
            %OPT_SYNTHESIS Construct an instance of this class
            %   Detailed explanation goes here

            obj@opt_manager_interface(sys);
                        
            if nargin == 1 || isempty(iqc_op)                
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
                nop = length(obj.sys.op);
                iqc_op = cell(nop);
                bind =obj.sys.bind;
                for i = 1:nop
                    nrep = sum(i==bind);
                    op_blank = obj.sys.op{i}.create_iqc_identity(nrep);
                    
                    iqc_op{i} = op_blank.factor();
                end
        end
        
        function obj = process_argument(obj,iqc_op)
            %PROCESS_ARGUMENT assign orders to the operators/IQCs
            
            %iqc_rob: IQCs representing the robust uncertainties

            if nargin > 1 && ~isempty(iqc_op)
                obj.iqc_op = iqc_op;
            else
                obj.iqc_op = obj.make_blank_iqc();
            end
        end

        function [diss] = index_specs(obj, alg_psi, iqc_op, specs)

            %INDEX_SPECS:  index into the performance specifications
            %
            %
            %   diss:   structure describing the problem
            %       plant:  system to control
            %       spec:   performance specification           
            %       target: whether the performance measure should be optimized
            %               true:  soft constraint (e.g. Schur complement
            %                                       formulation)
            %               false: hard constraint
            
            %TODO: maybe this should go inside the (system), not (manager)?
            
            if nargin < 4
                specs = obj.specs;
            end

            if nargin < 5
                target_ind = 0;
            end


            diss = cell(length(specs), 1);
            %determine the indices for each performance specification
            for i = 1:length(specs)
                
                      
                sp = specs{i};
                iwp_iqc = (1:(iqc_op.nw))';
                ir_iqc_first = (1:(iqc_op.np))';


                count_iqc_in = (iqc_op.nw);
                count_iqc_out = (iqc_op.np);

                if isempty(sp.izp) || isempty(sp.iwp)
                    ir_iqc_first_r =[];
                    iw_iqc_first_r = [];
                else
                    iw_iqc_first_r = count_iqc_in + (1:sp.iwp);
                    count_iqc_in = count_iqc_in + sp.iwp;

                    ir_iqc_first_r = count_iqc_out  + (1:sp.izp);
                    count_iqc_out = count_iqc_out + sp.izp;
                end

                iwp_iqc = [iwp_iqc; iw_iqc_first_r];

                ir_iqc0 = [ir_iqc_first; ir_iqc_first_r];
               
                sp_ind_w = iwp_iqc;
                sp_ind_r = ir_iqc0;


                
                    


                if iscell(alg_psi)
                    %TODO: change to genplant_poly type?
                    nwr = alg_psi{1}.nz;
                    nww = alg_psi{1}.nw;                    
                else
                    nwr = alg_psi.nz;
                    nww = alg_psi.nw;                    
                end

                nu = obj.sys.nu;
                ny = obj.sys.ny;
                E_r = blkdiag(full(sparse(1:length(sp_ind_r), sp_ind_r, ones(1, length(sp_ind_r)), length(sp_ind_r), nwr)), eye(ny));
                E_w = blkdiag(full(sparse(1:length(sp_ind_w), sp_ind_w, ones(1, length(sp_ind_w)), length(sp_ind_w), nww)), eye(nu));

                

                %enforce squareness in the performance specs?

                

                %nonminimal representation of the multiplier-extended plant

                %TODO: write fancier index code?
                

                n2 = alg_psi.dump_dim();
                n2.nwp = length(sp.iwp);
                n2.nzp = length(sp.izp);

                if iscell(alg_psi)

                    alg_screen = cell(size(alg_psi));
                    for j = 1:length(alg_screen)
                        alg_screen{j} = E_r * alg_psi{j} * E_w;
                    end

                    %TODO: write this part: cells/genplant poly
                else

                    alg_screen_P = E_r * alg_psi.ss * E_w;
                    alg_screen = genplant(alg_screen_P, n2);

                    
                    
                end


                diss{i} = struct('iqc_rob', iqc_op, ...
                    'spec', sp);
                diss{i}.plant = alg_screen;
                % %need to permute the entries of Mdiag for the partition





                %TODO: this may run into trouble if one entry has an X.
                %performance with dynamic multipliers?
            
                % diss{i} = struct('plant', alg_screen, 'M', M, 'X', iqc_op.X, ...
                    % 'spec', sp);
            end

        end

        %% extract the solution                   
        function  sol = process_recovery(obj, sol, lmi_out, alg_psi)
            %PROCESS_RECOVERY recover the controller from the solution
            
            sol = obj.lmi.process_recovery(sol, lmi_out, alg_psi);
            
            % iqc_rec = cell(size(obj.iqc_op));
            % for i = 1:length(obj.iqc_op)
            %     if isnumeric(obj.iqc_op{i})
            %         %the Same oracle (m=L, known linear transformation)
            %         iqc_rec{i} = obj.iqc_op{i};
            %     else
            %         iqc_rec{i} = obj.iqc_op{i}.recover(lmi_out);
            %     end
            % 
            % end
            % 
            % sol.iqc = iqc_rec;
        end
    end
end

