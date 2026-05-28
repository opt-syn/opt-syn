classdef lmi_synthesis_lti_reduced_order < lmi_synthesis_lti
    %LMI_SYNTHESIS_LTI synthesis LMIs for algorithmic interconnections
    %involving linear-time-invariant (LTI) networks and controllers
    %
    % w(k) \in F(z(k))
    %
    % [x(k+1)] = [A    B     Bp     Bu  ][x(k)]   state transition
    % [z(k)  ] = [Cz   Dzw   Dzwp   Dzu ][w(k)]   input to oracle
    % [zp(k) ] = [Czp  Dzpw  Dzpwp  Dzpu][wp(k)]  performance specification
    % [y(k)  ] = [Cy   Dyw   Dywp   Dyu ][u(k)]   input to controller
    %
    % performance specification: wp -> zp from (spec)
    %
    %   Implemented:
    %
    %   TODO:
    %       stability
    %       e2e
    %       quad
    %       p2p
    %       h2      
    %       e2p
    %       
    %
    
    methods
        function obj = lmi_synthesis_lti_reduced_order(sys, config)
            %LMI_SYNTHESIS_LTI Construct an instance of this class
            %   Detailed explanation goes here
            obj@lmi_synthesis_lti(sys, config);
        end       
        
        %% reduced-order control indexers
        function [vars, cons] = create_vars(obj, vars, cons, alg_psi, specs)
            %CREATE_VARS create the variables for the problem

            [vars, cons]  = create_vars@lmi_synthesis_interface(obj, vars, cons, alg_psi, specs);

           
            vars.diss.GS = obj.Pibar(vars.reg);
           
        end

        function Pb = Pibar(obj, vars_reg)
            %similarity transformation for optimization over Pi
            %used in regulator (reduced-order)            
            n = length(obj.sys.P.A);
            np = ssize(vars_reg.Pi, 1);
            ns = ssize(vars_reg.Pi, 2);
            
            Pb = [[eye(n-np), zeros(n-np, n+ns)];
                   [zeros(np, n-np), eye(np), vars_reg.Pi]];

        end

        function Ph = Pihat(obj, vars_reg)
            %similarity transformation for reduced-order control
            %used in regulator (reduced-order)
            n = length(obj.sys.A);
            np = ssize(vars_reg.Pi, 1);
            ns = ssize(vars_reg.Pi, 2);
            
            Ph = [Pb; [zeros(ns, n), eye(ns)]];
        end

        function cons = con_spread_single(obj, cons, GX, GY)
            %CON_SPREAD_SINGLE increase numerical conditioning by separating the 
            %primal and dual blocks
            % np = ssize(GX, 1);
            % spr = obj.config.tol.spread+1;
            % cons_PH = [GX, (spr)*eye(np); (spr)*eye(np), GY];
            % cons = append_lmi(cons, cons_PH, obj.LMILAB);

            cons = [];
        end
        
        function supp_b = supply_block(obj, sys_cl, quad)
            error('reduced_order: supply block not yet supported')
        end

        function stor_b = storage_block(obj, sys_cl, quad, G1, G2)
            error('reduced_order: storage block not yet supported')
        end

        function [sys_cl, U_cl, V_cl] = system_closed_loop(obj, P, vars_diss, vars_reg, vars_K);
            error('reduced_order: system closed loop block not yet supported')
        end
        
        %% recovery

        function [K_nofeed] = recover_subcontroller_warp(obj, P_trans, vars_rec)

            %RECOVER_SUBCONTROLLER_WARP recover the nonlinearly warped
            %controller 
            %dynamics and indexers


            %for debugging
            % G = obj.get_storage(sol.vars.diss, sol.vars.reg);

            %this is the (nonlinearly-warped) system that is certified as
            %possessing the desired performance and robustness
            %specifications
            % sys_cl = obj.system_closed_loop(P_trans, sol.vars.diss, sol.vars.reg, sol.vars.K);

            % sys_cal = ss(G \ Acl, G \ Bcl, Ccl, Dcl, 1);


            error('reduced_order: recovery not yet supported')
            [A, B, C, D] = ssdata(P_trans);

            iz = [P_trans.index_z(), P_trans.index_zp()];
            iw = [P_trans.index_w(), P_trans.index_wp()];
            iu = P_trans.index_u();
            iy = P_trans.index_y();           

            nz = length(iz);
            nw = length(iw);
            nu = length(iu);
            ny = length(iy);
                        

            [Ak, Bk, Ck, Dk] = obj.recover_K_from_elim(vars_rec);                        

            S = (vars_rec.diss.GS);
            n = ssize(Ak, 1);

            Y = vars_rec.diss.GY;
            X = vars_rec.diss.GX;


            J = S - X * Y;
            [Up, Sig, Vp] = svd(J);

            % U = Up*Sig;
            ssig = sqrt(Sig);
            srsig = diag(1./(diag(ssig)));


            V = Vp*ssig;
            U = Up*ssig;

            Uinv = srsig*Up';
            Vinv = srsig*Vp';


            %similarity transformation

            % 
            % %get right-side entries
            % I = eye(n);
            % Z1 = (Vinv*(I - X * Y')')';
            % Z2 = (Vinv* (-U * Y')')';
            % 
            % Z34 = [X, Z1; U, Z2] \ [zeros(n); eye(n)];
            % 
            % Z3 = Z34(1:n, :);
            % Z4 = Z34((n+1):end, :);
            % 
            % T = [eye(n), Y'; zeros(n), V'];
            % Ti = [eye(n), -Y' * Vinv'; zeros(n), Vinv'];
            % 
            % SimG = [X, Z1; U, Z2];
            % SimGi = [Y', Z3; V', Z4];

            %controller recovery

            Lblock = [Uinv, -Uinv*X*B(:, iu);
                zeros(nu, size(V, 2)), eye(nu)];

            LblockI = [U, X*B(:, iu); 
                zeros(nu, size(V, 2)), eye(nu)];

            Cblock = [Ak - X*A*Y, Bk;
                Ck, Dk];

            RblockI = [V' , zeros(size(V, 2), ny);
                C(iy, :)*Y, eye(ny)];

            Rblock = [Vinv', zeros(size(V, 2), ny);
                -C(iy, :)*Y*Vinv', eye(ny)];


            % Kblock0 = inv(LblockI)* (Cblock) * inv(RblockI);
            % Kblock1 = LblockI) \ Cblock) * inv(RblockI);
            % Kblockinv = RblockI') \ LblockI) \ Cblock)')';

            Kblock = Lblock * Cblock * Rblock;

            %extraction and exponential weighting
            Ac = Kblock(1:n, 1:n);
            Bc = Kblock(1:n, n+1:end);
            Cc = Kblock(n+1:end, 1:n);
            Dc = Kblock(n+1:end, n+1:end);

            K_nofeed_full = ss(Ac, Bc, Cc, Dc, 1);
            K_nofeed = K_nofeed_full;            
        end



        

    end

    
end


