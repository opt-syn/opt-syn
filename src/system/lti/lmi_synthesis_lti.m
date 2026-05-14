classdef lmi_synthesis_lti < lmi_synthesis_interface
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
        function obj = lmi_synthesis_lti(sys, config)
            %LMI_SYNTHESIS_LTI Construct an instance of this class
            %   Detailed explanation goes here
            obj@lmi_synthesis_interface(sys, config);
        end       
        
        
        %% definition of variables and helpers
        function [vars_diss, cons]= create_vars_storage(obj, cons, alg_psi, name)
            %create_vars_storage create variables for the dissipation
            %constraints
            %
            %Input:
            %   cons:       accumulated constraints
            %   alg_psi:    the filtered algorithmic interconnection
            %   name:       a name for the variable

            if nargin < 4
                name = [];
            end

            
            [GX, GY, cons] = obj.define_storage_G(cons, alg_psi, name);
            n = ssize(GX, 1);
            vars_diss= struct('GX', GX, 'GY', GY, 'S', eye(n));

        end

        %% stability (for testing)
        function [cons, objective, con_M] = stability_passive(obj, vars, cons, diss)
        % function [cons, objective, con_M] = stability(obj, vars, cons, diss)
            %certification of exponential stability
            
            G = obj.get_storage(vars.diss, vars.reg);

            %IMPORTANT!
            %hook up the internal model
            %(maybe it should happen at a higher level?)
            P = obj.reg.connect_model(diss.plant, diss.rho);

            sys_cl = obj.system_closed_loop(P, vars.diss, vars.reg, vars.K);

            %only do this if the system is passive

            n = ssize(sys_cl.A, 1);
            nw = ssize(sys_cl.B, 2);
            nz = ssize(sys_cl.C, 1);
            
            % nt = length(ind_pos);

            dyn_block =  [G,  sys_cl.A, sys_cl.B;
            sys_cl.A', G, zeros(n, nz);
            sys_cl.B', zeros(nz, n), zeros(nz)];

            %supply block
            %       sp = [zeros(n), zeros(n), zeros(n, nz);
            % zeros(n), zeros(n), zeros(n, nz); %vars.ga * eye(n)
            % zeros(nz, n), Ccl, Dcl + eye(nz)*p.opts.pass_tol] * (-0.5);

            dissI = eye(nz)*obj.config.tol.input_diss;
            % dissI = zeros(nz);
              sp = [zeros(n), zeros(n), zeros(n, nz);
            zeros(n), zeros(n), zeros(n, nz); %vars.ga * eye(n)
            zeros(nz, n), sys_cl.C, sys_cl.D + dissI] * (-0.5);
            
            cost_block = sp + sp';
    
            con_M = dyn_block + cost_block;
    
            objective = 0;            


            sM = ssize(con_M,1);
            cons = append_lmi(cons, con_M - obj.config.tol.M*eye(sM), obj.LMILAB); 

            %impose sign constraint
            %change this up
            % cons = obj.con_terminal(G, cons, [], diss.iqc_rob);
        end


        %% Quadratic performance (infinite horizon)        
        function [cons, objective, con_M] = quad(obj, vars, cons, diss)
            %QUAD: certificate of infinite-horizon quadratic performance
            
            %get the variables of the problem
            G = obj.get_storage(vars.diss, vars.reg);
            
            %IMPORTANT!
            %hook up the internal model
            %(maybe it should happen at a higher level?)
            P = obj.reg.connect_model(diss.plant, diss.rho);

            sys_cl = obj.system_closed_loop(P, vars.diss, vars.reg, vars.K);
            
            %index the quadratic specification
            vars_spec = vars.spec{diss.spec.id};
            M_quad = -obj.merge_spec_M(diss.iqc_rob, diss.spec, vars_spec);


            if isempty(diss.spec.izp)
                ind_p = 1:(diss.iqc_rob.nz);
                ind_q = diss.iqc_rob.nz + (1:(diss.iqc_rob.nw));
            else
                ind_p = 1:(diss.iqc_rob.nz + diss.spec.izp);
                ind_q = (diss.iqc_rob.nz + diss.spec.izp) + (1:(diss.iqc_rob.nw + diss.spec.iwp));
            end
            

            quad = obj.quad_objective(M_quad, ind_p, ind_q);
            
            
            %formulation from ParDynSyn notes (parametric dynamic
            %synthesis)


            %the quadratic objective
            supp_b = obj.supply_block(sys_cl, quad);

            %the storage
            stor_b = obj.storage_block(sys_cl, quad, G, G);

            %the dynamics
            dyn_b = obj.dynamics_block(sys_cl, quad);
            
            %wrap it all together
            objective = 0;

            con_M = stor_b + supp_b + dyn_b;


            sM = ssize(con_M,1);
            cons = append_lmi(cons, con_M - obj.config.tol.M*eye(sM), obj.LMILAB); 

            %impose sign constraint            
            cons = obj.con_terminal(G, cons, [], diss.iqc_rob);
        end        



        function [cons, objective, con_M] = e2e_target(obj, vars, cons, diss)
            %E2E_TARGET: use a Schur complement to minimize the energy to
            %energy gain of the transfer function


            
                       %get the variables of the problem
            G = obj.get_storage(vars.diss, vars.reg);
            
            %IMPORTANT!
            %hook up the internal model
            %(maybe it should happen at a higher level?)
            P = obj.reg.connect_model(diss.plant, diss.rho);

            sys_cl = obj.system_closed_loop(P, vars.diss, vars.reg, vars.K);
            
            %index the quadratic specification
            vars_spec = vars.spec{diss.spec.id};
            M_quad = -diss.iqc_rob.M;            


            ind_p = 1:(diss.iqc_rob.nz);
            ind_q = diss.iqc_rob.nz + (1:(diss.iqc_rob.nw));
            
            
            mu_l2 = vars.spec{diss.spec.id}.mu_l2;

            quad_rob = obj.quad_objective(M_quad, ind_p, ind_q);

            %adapt the quadratic objecitive for the e2e target
            nwp = length(diss.spec.iwp);
            nzp = length(diss.spec.izp);
            
            Q_e2e = -eye(nzp) * mu_l2;
            T_e2e = eye(nwp);
            S_e2e = zeros(nzp, nwp);
            U_e2e = -eye(nzp) * mu_l2;

            Q_new = blkdiag(quad_rob.Q, Q_e2e);
            T_new = blkdiag(quad_rob.T, T_e2e);
            S_new = blkdiag(quad_rob.S, S_e2e);
            U_new = blkdiag(quad_rob.U, U_e2e);

            quad = struct('Q', Q_new, 'T', T_new, 'S', S_new, 'U', U_new);

            
            %formulation from ParDynSyn notes (parametric dynamic
            %synthesis)


            %the quadratic objective
            supp_b = obj.supply_block(sys_cl, quad);

            %the storage
            stor_b = obj.storage_block(sys_cl, quad, G, G);

            %the dynamics
            dyn_b = obj.dynamics_block(sys_cl, quad);
            
            %wrap it all together
            objective = mu_l2;
            con_M = stor_b + supp_b + dyn_b;


            sM = ssize(con_M,1);
            cons = append_lmi(cons, con_M - obj.config.tol.M*eye(sM), obj.LMILAB); 

            %impose sign constraint            
            cons = obj.con_terminal(G, cons, [], diss.iqc_rob);
            
            
           
        end

        %% Peak-to-Peak norm (at each finite horizon)

        function [cons, objective, con_M] = p2p(obj, vars, cons, diss)
            %p2p: certificate of peak to peak induced norm
            %
            % sup norm(zp, 2) / norm(wp, 2) <= objective

            % verification by Theorem 4 of https://www.sciencedirect.com/science/article/pii/S2405896323008194

            %TODO: fix exponential rate here
            %storage matrix

            error('LTI synthesis: p2p target not yet supported')
            
        end


        

    end

    
end


