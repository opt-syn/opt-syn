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
        

        function [diss, cons] = build_plant(obj, cons)
            %BUILD_PLANT: form the plant to be used for analysis           

            %
            %Output:
            %   diss:   structure for the dissipation constraint
            %       plant: system 
            %       M:      running cost
            %       X:      terminal cost
            %       rho:    scaling/discount factor
            alg = obj.sys.get_alg;

            %sort based on the bind
            nop = length(obj.sys.bind);
            [~, perm] = sort(obj.sys.bind);
            
            c = obj.sys.P.nw/length(obj.sys.op);
            P = eye(nop);
            P(:, perm) = P;
            P = kron(P, eye(c));

            %TODO: allow for dynamical uncertainties in the w channel
            Pwp = blkdiag(P', eye(obj.sys.P.nwp));
            Pzp = blkdiag(P, eye(obj.sys.P.nzp));

            alg_perm = Pzp * alg * Pwp;

            %TODO: get rid of the same (m=L) oracles
            
            iqc_op = obj.iqc_op_all();

            alg_loop = lft(iqc_op.loop, alg_perm, nop, nop);


            if isempty(obj.specs)
                specs = opt_performance('stab', 1 - 1e-4);
            else
                specs = obj.specs;
            end

            
            diss = cell(length(specs), 1);
            for i = 1:length(specs)
                sp = specs{i};
                % iqc_spec = sp.iqc;
                
                [vars, cons, iqc_spec] = sp.create_iqc(cons);

                %selector matrices
                sp_w_ind = obj.sys.P.nw + sp.iwp;
                sp_z_ind = obj.sys.P.nz + sp.izp;
                
                E_wp = full(sparse(1:length(sp.iwp), sp.iwp, ones(1, length(sp.iwp)), length(sp.iwp), obj.sys.P.nwp));
                E_zp = full(sparse(1:length(sp.izp), sp.izp, ones(1, length(sp.izp)), length(sp.izp), obj.sys.P.nzp));

                E_w = blkdiag(eye(obj.sys.P.nw), E_wp);
                E_z = blkdiag(eye(obj.sys.P.nz), E_zp);
                
                % TODO: allow for performance IQCs with general loop
                % transformations, right now performance doesn't allow for
                % this. Only loop transformations on the oracles are
                % permitted.
                 
                alg_screen = E_w * alg_loop * E_z;

                iqc = blkdiag(iqc_op, iqc_spec);

                psi = iqc.get_psi();

                I = ss(eye(nop + length(sp.iwp)));

                %VERY IMPORTANT: [G; I], not blkdiag(G, I) (like I
                %previously had, embarassing bug was here, 21.04.2026)
                GI = [alg_screen; I];
                alg_psi = psi * GI;
    
                diss{i} = struct('plant', alg_psi, 'M', iqc.M, 'X', iqc.X, 'rho', sp.rho, 'np', iqc.np, 'vars', vars);

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
            %BUILD_PROGRAM set up the analysis problem
            
            
            cons = obj.cons;

          
            [diss, cons] = obj.build_plant(cons);
            ndiss = length(diss);
            
            [vars_diss, cons] = obj.create_vars_diss(diss, cons);

            cons = obj.coeff_normalize(cons);

           
           
            
            %incorporate parameters afterwards
            param = [];
            for i = 1:ndiss
                [con_M, con_X] = build_dissipation(obj, vars_diss, diss{i}, param);

                
                if obj.opts.impose_X && (i==1)
                    sx = ssize(con_X, 1);
                    cons = append_lmi(cons, con_X - eye(sx)*obj.tol.X, obj.LMILAB);
                end
                sm = ssize(con_M, 1);
                cons = append_lmi(cons, con_M - eye(sm)*obj.tol.M, obj.LMILAB);
            end

            vars = obj.vars;
            vars.diss = vars_diss;

            vars.specs = cell(1, length(obj.specs));
            for i = 1:length(obj.specs)
                vars.specs{i} = diss{i}.vars;
            end
        end

        function [vars_diss, cons]= create_vars_diss(obj, diss, cons)
            %CREATE_VARS_DISS create variables for the dissipation
            %constraints

            n = length(diss{1}.plant.A);
            G = lmim('G', n, n, 'sym');
            %TODO: TESTING
         %    G = [1.0000         0
         % 0    6.0500];

            % ga = lmim('ga', 1, 1);

            vars_diss= struct('G', G);
            % vars_diss= struct('G', G, 'ga', ga);
            
            % cons  = append_lmi(cons, ga, obj.LMILAB);


                        %norm bound 
                                    % cons = append_lmi(cons, ga - trace(G), obj.LMILAB);
% 
            if obj.tol.G_max < Inf    
                %issue in the bounding?
                % cons = append_lmi(cons, obj.tol.G_max*eye(2*n) + [zeros(n), G; G', zeros(n)], obj.LMILAB);

                % cons = append_lmi(cons, obj.tol.G_max*eye(n)  - G, obj.LMILAB);
                % cons = append_lmi(cons, obj.tol.G_max*eye(n)  + G, obj.LMILAB);

                cons = append_lmi(cons, obj.tol.G_max*eye(n)  - G, obj.LMILAB);
                cons = append_lmi(cons, obj.tol.G_max*eye(n)  + G, obj.LMILAB);

            end

            % cons = append_lmi(cons, G - obj.tol.G*eye(n), obj.LMILAB);
        

            
            
        end


        function  sol = process_recovery(obj, sol, lmi_out);
            %PROCESS_RECOVERY post-process the solution
            
            %TODO: fill this in

            iqc_rec = cell(size(obj.iqc_op));
            for i = 1:length(obj.iqc_op)
                iqc_rec{i} = obj.iqc_op{i}.recover(lmi_out);

            end

            sol.iqc = iqc_rec;
            
            




        end
        
        
        
        function [con_M, con_X] = build_dissipation(obj, vars_diss, diss, param)
            %FORM_DISSIPATION: the dissipation relation for IQC synthesis

            %TODO: special calls for the h infinity and h2 structures
            %for convexification

            G = vars_diss.G;

            %storage routines
            % 
            G0 = [1.0000         0
            0    6.0500];
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

