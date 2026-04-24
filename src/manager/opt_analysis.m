classdef opt_analysis < opt_manager_interface
    %OPT_ANALYSIS  analysis of optimization algorithms
    %
    % iterative procedure to find a point beta satisfying
    % the fixed-point equation 
    %               0 \in sum_i F_i(\beta).
    %
       
    methods
        function obj = opt_analysis(sys)
            %OPT_ANALYSIS Construct an instance of this class            
           
            obj@opt_manager_interface(sys);

        end
        
        %% define IQCs for the operators
        function [obj, vars, cons] = oracle_order(obj,order, ind)
            %ORACLE_ORDER: set the orders of the IQCs
            nop = length(obj.sys.op);
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

        function [iqc, m_same, ind_same] = iqc_op_all(obj)
            %IQC_OP_ALL: all iqcs for the operators
            iqc = {};
            m_same = [];
            ind_same = [];
            same_count = 0;

            for i = 1:length(obj.iqc_op)
                if ~obj.sys.op{i}.same
                    %block diagonal of the iqc
                    if isempty(iqc)
                        iqc = obj.iqc_op{i};
                    else
                        iqc = blkdiag(iqc, obj.iqc_op{i});
                    end
                    same_count = same_count + obj.iqc_op{i}.nw;
                    
                else
                    %treat the m=L case separately
                    m_same = blkdiag(m_same, obj.iqc_op{i});

                    ind_same = [ind_same, same_count + (1:length(m_same))];
                    same_count = same_count + length(m_same);
                end            
            end

        end


        %% build the system
        function [alg_psi, iqc_op, alg_loop] = build_plant(obj, cons)
            %BUILD_PLANT: form the plant to be used for analysis           

            %
            %Output:
            %   alg_psi:    plant with filters (psi)
            %   alg_loop:   plant without filters, but after loop
            %               transformation (should be stable)
            %   iqc_op:     iqcs for the robust uncertainties
            alg = obj.sys.get_alg;

            %sort oracles based on the bind (exposure of repeated
            %nonlinearities)
            nop = length(obj.sys.bind);
            [~, perm] = sort(obj.sys.bind);
            
            c = obj.sys.op{1}.c;
            P = eye(nop);
            P(:, perm) = P;
            P = kron(P, eye(c));
            
            Pwp = blkdiag(P', eye(obj.sys.P.nwp));
            Pzp = blkdiag(P, eye(obj.sys.P.nzp));

            alg_perm = Pzp * alg * Pwp;

            
            %identify and get rid of the same (m=L) oracles   

            %use an explicit substitution

            [iqc_op, m_same, ind_same] = obj.iqc_op_all();

            ind_diff = setdiff(1:obj.sys.P.nw, ind_same);
            Pd = eye(nop);
            Pd(:, [ind_same, ind_diff]) = Pd;
            n_same = length(ind_same);

            Pwp2 = blkdiag(Pd', eye(obj.sys.P.nwp));
            Pzp2 = blkdiag(Pd, eye(obj.sys.P.nzp));
            alg_perm_same = Pzp2 * alg_perm * Pwp2;

            alg_perm_m = lft(m_same, alg_perm_same, n_same, n_same);
            
            %no loop transformations in performance
            loop = iqc_op.loop;
            nloop = length(loop)/2;
            alg_loop = lft(loop, alg_perm_m, nloop, nloop);

            %form the system
            I = ss(eye(size(alg_loop.D, 2)));
            GI = [alg_loop; I];

            obj.sys.P.nwp;
            
            Psi1 = iqc_op.Psi1;
            Psi2 = iqc_op.Psi2;
            I_zp = eye(obj.sys.P.nzp);
            I_wp = eye(obj.sys.P.nwp);

            psi = blkdiag(Psi1, I_zp, Psi2, I_wp);
           

            alg_psi = psi * GI;           
          
        end

        function [diss] = index_specs(obj, alg_psi, iqc_op, specs)

            %INDEX_SPECS:  index into the performance specifications
            %
            %
            %now index alg_psi into its performance specifications
            if nargin < 4
                specs = obj.specs;
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

                ir_iqc = [ir_iqc_first; ir_iqc_first_r];
                ir_iqc = [ir_iqc; ir_iqc + (iqc_op.np + obj.sys.P.nwp )];

                sp_ind_w = iwp_iqc;
                sp_ind_r = ir_iqc;
                % nww = length(sp_ind_w);
                % nwr = length(sp_ind_r);
                [nwr, nww] = ssize(alg_psi.D);
                %enforce squareness in the performance specs?
                E_w = full(sparse(1:length(sp_ind_w), sp_ind_w, ones(1, length(sp_ind_w)), length(sp_ind_w), nww));
                E_r = full(sparse(1:length(sp_ind_r), sp_ind_r, ones(1, length(sp_ind_r)), length(sp_ind_r), nwr));


                %nonminimal representation
                alg_screen = E_r * alg_psi * E_w;


                %need to permute the entries of Mdiag for the partition
                n1 = iqc_op.np;
                m1 = iqc_op.nq;
                n2 = length(sp.iwp);
                m2 = length(sp.izp);
                [M] = outer_blkdiag(iqc_op.M, sp.supply, n1, m1, n2, m2);
                % Mdiag = blkdiag(iqc_op.M, sp.supply);




                %TODO: this may run into trouble if one entry has an X.
                %performance with dynamic multipliers?
            
                diss{i} = struct('plant', alg_screen, 'M', M, 'X', iqc_op.X, ...
                    'spec', sp);
            end

        end

        function con_X = con_terminal(obj, vars, iqc_op)
            
            %terminal cost constraint (nonnegativity on G)
            X = iqc_op.X;
            G = vars.diss.G;

            nf = ssize(X);
            n = ssize(G, 1);
            Ef = [eye(nf); zeros(n-nf, nf)];

            X_f = Ef * X * Ef';
            con_X = G + X_f;
            
        end
        
        
        function [con_M] = con_dissipation(obj, vars, diss, param)
            %CON_DISSIPATION: the dissipation relation for IQC synthesis

            %TODO: special calls for the h infinity and h2 structures
            %for convexification

            if nargin < 4
                param = [];
            end

            G = vars.diss.G;

            %storage routines
            G_current = G;
            G_next = G;

            Gblock = blkdiag(diss.spec.rho^2 * G_current, -G_next);

            [n, m] = size(diss.plant.B);  


                Ablock = [eye(n), zeros(n, m);
                    diss.plant.A, diss.plant.B];
    
                Cblock = [diss.plant.C, diss.plant.D];
               
                                
                Center_block = blkdiag(Gblock, -diss.M);
                Outer_block = [Ablock; Cblock];

                con_M = (Outer_block' * Center_block * Outer_block);
            

                sys_block = Ablock' * Gblock * Ablock;
                
                supp_block = -Cblock' * diss.M * Cblock;
        end



        function objective = get_objective(obj, vars)
            %GET_OBJECTIVE determine the objective for optimization
            %TODO: fill in specification

            % objective = vars.diss.ga;
            % objective = trace(vars.diss.G);
            objective = 0;
            % objective = 0;

            % for i = 1:length(obj.specs)
            %     if strcmp(obj.specs{i}.type, 'finite_l2')
            %         objective = vars.specs{i}.mu_l2;
            %     end
            % end
        end

        function [vars, cons] = build_program(obj, specs)
            %BUILD_PROGRAM set up the algorithm analysis problem
            
            if nargin < 2
                specs = obj.specs;
            end
            
            cons = obj.cons;
            vars = obj.vars;
          
            % [diss, cons] = obj.build_plant(cons);
            %alg_loop: used for debugging. The algorithm after signal 
            % transformationsbefore, but before cascade by the filters    

            [alg_psi, iqc_op, alg_loop] = obj.build_plant(cons);
            [vars.diss, cons] = obj.create_vars_storage(alg_psi, cons);
            [vars.spec, cons] = obj.create_vars_spec(specs, cons);

            
            %the dissipation can change
            [vars, cons] = build_dissipation(obj, vars, cons, alg_psi, iqc_op, specs);

        end

        function [vars, cons] = build_dissipation(obj, vars, cons, alg_psi, iqc_op, specs)
            %BUILD_DISSIPATION: form the dissipation relations for the
            %system (at the current set of specifications)
           
            [diss] = obj.index_specs(alg_psi, iqc_op, specs);
            ndiss = length(diss);

            %dissipation relations
            for i = 1:ndiss
                con_M = con_dissipation(obj, vars, diss{i});                                 
                    
                                    
                sm = ssize(con_M, 1);
                cons = append_lmi(cons, con_M - eye(sm)*obj.tol.M, obj.LMILAB);
            end

            %terminal cost/sign constraints
            if obj.opts.impose_X
                %for infinite-horizon performance measures (l2 norm, h2
                    %norm), terminal costs and sign constraints on the
                    %storage function are not required. For finite-horizon
                    %specifications (e.g. invariance, peak-to-peak), they
                    %are needed.

                
                con_X = obj.con_terminal(vars, iqc_op);
                sx = ssize(con_X, 1);
                cons = append_lmi(cons, con_X - eye(sx)*obj.tol.X, obj.LMILAB);
            end

        end

        function [vars_diss, cons]= create_vars_storage(obj, alg_psi, cons)
            %create_vars_storage create variables for the dissipation
            %constraints

            n = length(alg_psi.A);
            G = lmim('G', n, n, 'sym');

            %ga
            % ga = lmim('ga', 1, 1, 'sym');

            vars_diss= struct('G', G);
            % vars_diss= struct('G', G, 'ga', ga);

            

            if obj.tol.G_max < Inf    
                %issue in the bounding?
                cons = append_lmi(cons, obj.tol.G_max*eye(n)  - G, obj.LMILAB);
                cons = append_lmi(cons, obj.tol.G_max*eye(n)  + G, obj.LMILAB);

                %lmim complains that the dimensions are wrong here.
                % cons = append_lmi(cons, ga*eye(n)  - G, obj.LMILAB);
                % cons = append_lmi(cons, ga*eye(n)  + G, obj.LMILAB);
            end


            
            
        end

        function  sol = process_recovery(obj, sol, lmi_out);
            %PROCESS_RECOVERY post-process the solution
            
            %TODO: fill this in

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

