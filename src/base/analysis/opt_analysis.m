classdef opt_analysis < opt_manager_interface
    %OPT_ANALYSIS Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        function obj = opt_analysis(sys)
            %OPT_ANALYSIS Construct an instance of this class
            %   Detailed explanation goes here
           
            obj@opt_manager_interface(sys);
        end
        
        function obj = oracle_order(obj,order, ind)
            %ORACLE_ORDER: set the orders of the IQCs
            nop = length(obj.sys.op);
            if nargin < 3
                ind = 1:nop;
            end
            
            for i = 1:nop
                if ismember(i, ind)
                    rep_curr = nnz(obj.sys.bind==i);
                    [iqc_curr, vars_curr,cons_curr] = obj.sys.op{i}.create_iqc(obj.cons, order{i}, rep_curr);
    
                    obj.iqc_op{i} = iqc_curr;
                    obj.vars.op{i} = vars_curr;
                    obj.cons = cons_curr;
                end
            end

            %normalize the coefficients for the filters
            obj.cons = obj.coeff_normalize(obj.cons);
        end


        function cons = coeff_normalize(obj, cons)
            %COEFF_NORMALIZE add constraint to normalize the psi multipliers
            
            nop = length(obj.sys.op);
            cs = 0;
            for i = 1:nop
                cs_curr = obj.sys.op{i}.csum_psi(obj.vars.op{i});
                cs = cs + cs_curr;

            end

            cons = append_lmi(cons, cs - nop*0.9, obj.LMILAB);
            cons = append_lmi(cons, -cs + nop*1.1, obj.LMILAB);


        end



        function [alg_psi, alg_all, data_spec] = build_plant(obj, cons)
            %BUILD_PLANT: form the plant to be used for analysis           

            %
            %Output:
            %   diss:   structure for the dissipation constraint
            %       plant: system 
            %       M:      running cost
            %       X:      terminal cost
            %       rho:    scaling/discount factor
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
            


            %we now have the loop terms for the operators. now go through 
            %the loop terms for the performance specifications


            if isempty(obj.specs)
                specs = {opt_performance('stab', 1 - 1e-4)};
            else
                specs = obj.specs;
            end

            iqc_all = iqc_op;            
            
            count_iqc_out = 0;
            data_spec = cell(length(specs), 1);
            
            for i = 1:length(specs)
                sp = specs{i};
                % iqc_spec = sp.iqc;
                
                [vars, cons, iqc_spec] = sp.create_iqc(cons);
                iqc_all = blkdiag(iqc_all, iqc_spec);
                iqc_curr = blkdiag(iqc_op, iqc_spec);


                if isempty(iqc_spec)
                    ind_iqc_out =[];
                else
                    ind_iqc_out = count_iqc_out  + 1:(iqc_spec.np );
                    count_iqc_out = count_iqc_out  + (iqc_spec.np);
                end
                data_spec{i} = struct('ind_spec', ind_iqc_out, 'iqc', iqc_spec, ...
                    'vars', vars, 'rho', sp.rho, 'n_same', n_same, 'iqc_with_op', iqc_curr);
            end

            %generate the entire system
            loop = iqc_all.loop;
            nloop = length(loop)/2;
            psi = iqc_all.get_psi();
            
            alg_all = lft(loop, alg_perm_m, nloop, nloop);

            
            I = ss(eye(nloop));

            %VERY IMPORTANT: [G; I], not blkdiag(G, I) (like I
            %previously had, embarassing bug was here, 21.04.2026)
            GI = [alg_all; I];
            alg_psi = psi * GI;           
        end

        function [diss] = index_specs(obj, alg_psi, data_spec )

            if isempty(obj.specs)
                specs = {opt_performance('stab', 1 - 1e-4)};
            else
                specs = obj.specs;
            end

            %now index alg_psi into its performance specifications
            diss = cell(length(specs), 1);
            for i = 1:length(specs)
                
                % iqc_spec = sp.iqc;
                n_same = data_spec{i}.n_same;

                sp_ind = data_spec{i}.ind_spec;
                nwp = length(sp_ind);
                %enforce squareness in the performanc specs?
                E_wp = full(sparse(1:length(sp_ind), sp_ind, ones(1, length(sp_ind)), length(sp_ind), nwp));
                E_zp = full(sparse(1:length(sp_ind), sp_ind, ones(1, length(sp_ind)), length(sp_ind), nwp));

                E_w = blkdiag(eye(obj.sys.P.nw - n_same), E_wp, eye(obj.sys.P.nw - n_same), E_wp);
                E_z = blkdiag(eye(obj.sys.P.nz- n_same), E_zp);
                
                 
                %nonminimal representation
                alg_screen = E_w * alg_psi * E_z;

                iqc = data_spec{i}.iqc_with_op;

                %TODO: this may run into trouble if one entry has an X.
                %performance with dynamic multipliers?
            
                diss{i} = struct('plant', alg_screen, 'M', iqc.M, 'X', iqc.X, ...
                    'rho', data_spec{i}.rho, 'np', iqc.np, 'vars', data_spec{i}.vars);
            end

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

        function objective = get_objective(obj, vars)
            %GET_OBJECTIVE determine the objective for optimization
            %TODO: fill in specification

            % objective = vars.diss.ga;
            % objective = -trace(vars.diss.G);
            % objective = 0;
            objective = 0;
            % for i = 1:length(obj.specs)
            %     if strcmp(obj.specs{i}.type, 'finite_l2')
            %         objective = -vars.specs{i}.mu_l2;
            %     end
            % end
        end

        function [vars, cons] = build_program(obj)
            %BUILD_PROGRAM set up the algorithm analysis problem
            
            
            cons = obj.cons;

          
            % [diss, cons] = obj.build_plant(cons);
            %alg_all: used for debugging. The algorithm after signal 
            % transformationsbefore, but before cascade by the filters            
            [alg_psi, alg_all, data_spec] = obj.build_plant(cons);
            [vars.diss, cons] = obj.create_vars_storage(alg_psi, cons);
            

            
            [vars, cons] = build_dissipation(obj, vars, cons, alg_psi, data_spec);

        end

        function [vars, cons] = build_dissipation(obj, vars, cons, alg_psi, data_spec)
            %BUILD_DISSIPATION: form the dissipation relations for the
            %system (at the current set of specifications)
           
            [diss] = obj.index_specs(alg_psi, data_spec);
            ndiss = length(diss);

            for i = 1:ndiss
                [con_M, con_X] = con_dissipation(obj, vars, diss{i});

                
                if obj.opts.impose_X && (i==1)
                    %for infinite-horizon performance measures (l2 norm, h2
                    %norm), terminal costs and sign constraints on the
                    %storage function are not required. For finite-horizon
                    %specifications (e.g. invariance, peak-to-peak), they
                    %are needed.
                    sx = ssize(con_X, 1);
                    cons = append_lmi(cons, con_X - eye(sx)*obj.tol.X, obj.LMILAB);
                end
                sm = ssize(con_M, 1);
                cons = append_lmi(cons, con_M - eye(sm)*obj.tol.M, obj.LMILAB);
            end



            vars.specs = cell(1, length(obj.specs));
            for i = 1:length(obj.specs)
                vars.specs{i} = diss{i}.vars;
            end

        end

        function [vars_diss, cons]= create_vars_storage(obj, alg_psi, cons)
            %create_vars_storage create variables for the dissipation
            %constraints

            n = length(alg_psi.A);
            G = lmim('G', n, n, 'sym');

            vars_diss= struct('G', G);

            if obj.tol.G_max < Inf    
                %issue in the bounding?
                cons = append_lmi(cons, obj.tol.G_max*eye(n)  - G, obj.LMILAB);
                cons = append_lmi(cons, obj.tol.G_max*eye(n)  + G, obj.LMILAB);

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
        
        
        
        function [con_M, con_X] = con_dissipation(obj, vars, diss, param)
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

            Gblock = blkdiag(diss.rho^2 * G_current, -G_next);


            [n, m] = size(diss.plant.B);  

            % %The true routines
            % STAB_TEST = false;

            % if  STAB_TEST
                
                % 
                % Ablock = [eye(n);
                % diss.plant.A];
                % Cblock = [diss.plant.C, diss.plant.D];
                % 
                % Center_block = Gblock;          
                % Outer_block = [Ablock];
                % 
                % 
                % % con_M = 10*eye(n) - G;
                % % con_M = eye(n);
                % % con_M = -G + 50*eye(n);
                % 
                % sys_block = Ablock' * Gblock * Ablock;
                
                % supp_block = -Cblock' * diss.M * Cblock;
            % else
                Ablock = [eye(n), zeros(n, m);
                    diss.plant.A, diss.plant.B];
    
                Cblock = [diss.plant.C, diss.plant.D];
               
                Center_block = blkdiag(Gblock, -diss.M);
                Outer_block = [Ablock; Cblock];

                con_M = (Outer_block' * Center_block * Outer_block);
            
            % end
            
            
            % con_M = diss.rho^2 * G_current - diss.plant.A' * G_next * diss.plant.A;

            %terminal cost constraint
            if isnumeric(diss.X)
                nf = length(diss.X);
            else
                nf = dim(diss.X, 1);
            end
            Ef = [eye(nf); zeros(n-nf, nf)];

            X_f = Ef * diss.X * Ef';
            con_X = G + X_f;
            
        end
    
    end
end

