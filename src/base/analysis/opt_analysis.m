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

        function [diss] = build_plant(obj)
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
            
            %TODO: make sure that the permutation is correct
            c = obj.sys.P.nw/length(obj.sys.op);
            P = eye(nop);
            P(:, perm) = P;
            P = kron(P, eye(c));

            Pwp = blkdiag(P', eye(obj.sys.P.nwp));
            Pzp = blkdiag(P, eye(obj.sys.P.nzp));

            alg_perm = Pzp * alg * Pwp; 
            
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
                iqc_spec = sp.iqc;

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

                alg_psi = psi * blkdiag(alg_screen, eye(nop + length(sp.iwp)));
    
                diss{i} = struct('plant', alg_psi, 'M', iqc.M, 'X', iqc.X, 'rho', sp.rho, 'np', iqc.np);

            end

            

        end

        function iqc = iqc_op_all(obj)
            %IQC_OP_ALL: all iqcs for the operators
            iqc = obj.iqc_op{1};
            for i = 2:length(obj.iqc_op)
                iqc = blkdiag(iqc, obj.iqc_op{i});
            end

        end


        function [vars, cons] = build_program(obj)
            %BUILD_PROGRAM set up the analysis problem
            diss = obj.build_plant();
            ndiss = length(diss);
            cons = obj.cons;

            n = length(diss{1}.plant.A);
            G = lmim('G', n, n, 'sym');

            cons = append_lmi(cons, G - obj.tol.G*eye(n), obj.LMILAB);
            for i = 1:ndiss
                [con_M, con_X] = build_dissipation(obj, G, diss{i});

                cons = append_lmi(cons, con_M, obj.LMILAB);
                cons = append_lmi(cons, con_X, obj.LMILAB);
            end

            vars = obj.vars;
            vars.diss = {'G', G};
        end

        function [con_M, con_X] = build_dissipation(obj, G, diss)
            %FORM_DISSIPATION: the dissipation relation for IQC synthesis

            %TODO: special calls for the h infinity and h2 structures
            %for convexification

            [n, m] = size(diss.plant.B);
            p = dim(diss.plant.C, 1);
                   

            Ablock = [eye(n), zeros(n, m);
                diss.plant.A, diss.plant.B];

            Gblock = blkdiag(diss.rho^2 * G, -G);

            store_block = Ablock' * Gblock * Ablock;
            % Cblock = [diss.plant.C, diss.plant.D;
                      % zeros(diss.np, n + m-diss.np), eye(diss.np)];
            Cblock = [diss.plant.C, diss.plant.D];

            supply_block = Cblock' * diss.M * Cblock;
            
            con_M = store_block - supply_block;

            %terminal cost constraint
            if isnumeric(diss.X)
                nf = length(diss.X);
            else
                nf = dim(diss.X, 1);
            end
            Ef = [eye(nf); zeros(n-nf, nf)];

            con_X = Ef' * G * Ef - diss.X;
            
        end
    
    end
end

