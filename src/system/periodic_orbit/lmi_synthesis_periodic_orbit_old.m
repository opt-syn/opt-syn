classdef lmi_synthesis_periodic_orbit < lmi_synthesis_lti
    %LMI_SYNTHESIS_PERIODIC_ORBIT synthesisLMIs for algorithmic interconnections
    %involving periodic linear networks and controllers
    %
    % Orbit structure on the periodicity
    %
    % w(k) \in F(z(k))
    %
    % [x(k+1)] = [A(k)    Bw(k)    Bwp(k)   Bu(k)  ][x(k)]   state transition
    % [z(k)  ] = [Cz(k)   Dzw(k)   Dzwp(k)  Dzu(k) ][w(k)]   output to oracle
    % [zp(k) ] = [Czp(k)  Dzpw(k)  Dzpwp(k) Dzpu(k)][wp(k)]  output to performance
    % [zp(k) ] = [Cy(k)   Dyw(k)   Dywp(k)  Dyu(k) ][u(k)]   output to controller
    %
    %A(k) = A(k+T) for some known time T
    %
    %instances of these algorithms include cyclic coordinate descent
    %methods. Periodic systems can also be unrolled into an LTI system
    %(monodromy methods): a single large LMI system rather than multiple 
    % coupled smaller LMI systems


    methods
        function obj = lmi_synthesis_periodic_orbit(sys,config)
            %LMI_SYNTHESIS_PERIODIC undefined
            %   undefined
            obj@lmi_synthesis_lti(sys, config);
        end

        %% definition of variables and helpers

        function ns = Nss(obj)
            %NSS: Number of subsystems            
            ns = obj.sys.Nss;
        end

        %% main call
        function [vars, cons, objective, con_M] = cons_dynamic(obj, vars, cons, diss)
            %CONS form the dissipation and sign constraints
            %
            %Input:
            %   vars:   variables of the problem        
            %   cons:   accumulated constraints
            %   diss:   structure describing the problem
            %       plant:  system to control
            %       spec:   performance specification           
            %       target: whether the performance measure should be optimized
            %               true:  soft constraint (e.g. Schur complement
            %                                       formulation)
            %               false: hard constraint            
            %
            %Output:
            %   cons:   accumulated constraints
            %   objective:  term to be minimized            


            %multiply the A and B entries of the plant by R
            %to accomplish the periodic orbit constructions


            n = size(diss.plant.A, 1);
            c = size(obj.sys.R, 1);

            diss.R = obj.sys.R;

            

            %Upper-levels: iterate over the systems
            [cons, objective, con_M] = obj.con_dynamic_single(vars, cons, diss);



        end

        function sol = process_recovery(obj, sol, lmi_out, alg_psi, diss)
            %recover the controller

            if nargin < 5
                diss = [];
            end

            %get the system with the internal model
            dissend = struct;
            dissend.plant = alg_psi;
            dissend.rho = sol.rho;
            dissend.R = obj.sys.R;
            P_trans =  obj.connect_model(dissend);

            %evaluate the variables
            [sol] = obj.recover_subcontroller(alg_psi, P_trans, sol);

            % sol.G = obj.get_storage(sol.vars.diss, sol.vars.reg);
        end

        function [sol] = recover_subcontroller(obj, alg_psi, P_trans, sol)
            %RECOVER_SUBCONTROLLER recover the subcontroller of the current
            %mode/control
            %
            %
            %Input:
            %
            %Output:
            %   K_feed: the subcontroller with direct feedthrough, before
            %           exponential discounting    
            %(not yet exponentially undiscounted, this happens later)

            vars_rec = sol.vars;
            rho = sol.rho;

            [K_nofeed, Gcal, Ycal] = recover_subcontroller_warp(obj, P_trans, vars_rec);

            model = obj.reg.get_model(1, vars_rec.reg);

            K_report = obj.K_alg_report(P_trans, K_nofeed, model, rho);

            sol.cert.alg_trans = K_report.alg_trans;
            sol.cert.alg = lft(obj.sys.P, K_report.K);
            sol.cert.model = K_report.model;           
            sol.K= K_report.K;
            sol.cert.K_sub = K_report.K_sub;
            sol.cert.Gcl = Gcal;
            sol.cert.Ycl = Ycal;

            sol.gain = obj.validate_recovery_gain(sol.cert.alg_trans, sol.cert.iqc_op_all);


        end

        function P_model = connect_model(obj, diss)
            reg_ind = 1;
            % reg_ind = obj.sys.Nss;
            P_model = obj.reg.connect_model(diss.plant, reg_ind, diss.rho);

            n = size(P_model.A);           
            c = size(diss.R, 1);
            Rkron = kron(eye(n/c), diss.R);

            p2 = P_model.P;
            p2.A = Rkron * p2.A;
            p2.B = Rkron * p2.B;

            P_model.P = p2;


        end

        %% Recovery
        % Call the recovery method to finalize the synthesis process
        % [cons, objective, con_M] = obj.recovery(vars, cons, diss);

        % function gain = validate_recovery_gain(obj, alg_trans, iqc_op_all)
        %     %VALIDATE_RECOVERY validate that the system obeys the stability
        %     %constraint (TODO: performance specs)
        % 
        %     %use the monodromy system to get specs
        %     n = alg_trans.dump_dim();
        % 
        %     sys_trans = obj.sys;
        %     sys_trans.K = [];
        %     sys_trans.P = alg_trans;
        %     sys_trans.P.P;
        %     %TODO: fix this
        %     alg_trans_lti = periodic_lift(sys_trans);
        % 
        %     P = alg_trans_lti.P; 
        % 
        %     c = obj.sys.op{1}.c;
        % 
        %     M = kron(iqc_op_all.iqc.M, eye(c));
        %     M = (M + M')/2;
        %     nw = floor(size(M, 1)/2);
        % 
        %     M11 = M(1:nw, 1:nw);
        %     M12 = M(nw + (1:nw), 1:nw);
        %     M22 = M(nw + (1:nw), nw + (1:nw));
        %     %is the constraint passive?
        %     is_passive = (norm(M11) + norm(M22) + norm(M12 - eye(nw)))==0;
        %     is_hinf = (norm(M11-eye(nw)) + norm(M22+eye(nw)) + norm(M12))==0;
        % 
        % 
        %     if is_passive
        %         gain_passive = -getPassiveIndex(-P, 'input');
        % 
        %         E=eye(nw);
        %         Tinf=[E sqrt(2)*E;sqrt(2)*E E];
        %         P_inf = lft(Tinf,P,nw,nw);
        % 
        %         gain_inf = norm(P_inf, 'inf');
        %     elseif is_hinf
        %         gain_inf = norm(P, 'inf');
        % 
        %         E=eye(nw);
        %         Tpass = [-E sqrt(2)*E;sqrt(2)*E -E];
        %         Ppass = lft(Tpass,P,nw,nw);
        % 
        %         gain_passive = -getPassiveIndex(-Ppass, 'input');
        %     else
        %         %TODO: advanced validation
        %         warning('Customized validation is not yet implemented')
        %         gain_inf = 0;
        %         gain_passive = 0;
        %     end
        % 
        %     gain = [gain_passive, gain_inf]; 
        % 
        % 
        %     % gain = validate_recovery_gain@lmi_synthesis_periodic(obj, alg_trans, iqc_op_all);
        % end
        % 
        % 

        function K_report = K_alg_report(obj, P_trans, K_nofeed, model, rho)
            %K_ALG_REPORT recover the algorithmic interconnection and the
            %controller
            
            D = P_trans.D;

            iz = [P_trans.index_z(), P_trans.index_zp()];
            iw = [P_trans.index_w(), P_trans.index_wp()];
            iu = P_trans.index_u();
            iy = P_trans.index_y();           

            nz = length(iz);
            nw = length(iw);
            nu = length(iu);
            ny = length(iy);

            %add the proper term by LFT
            D22 = D(iy, iu);
            Dfeed = zeros(nz+ny, nw+nu);            
            Dfeed(nz+1:end, nw+1:end) = D22;



            P_trans_nofeed = P_trans.ss;
            P_trans_nofeed.D = P_trans_nofeed.D - Dfeed;

            T_feed = [zeros(nu, ny), eye(nu); eye(ny), -D22];

            K_feed = lft(T_feed, K_nofeed, nu, ny);
            % K_feed_full = lft(T_feed, K_nofeed_full, nu, ny);

            
            K_feed = obj.name_K_feed(K_feed);



            alg_trans = lft(P_trans, K_feed);
            alg_trans_nofeed = lft(P_trans_nofeed, K_nofeed);


            K_sub= rhotrafo(K_feed, 1/rho);



            %account for the periodicity
            n = size(K_sub.A, 1);           
            c = size(obj.sys.R, 1);
            Rkron = kron(eye(n/c), obj.sys.R);

            p2 = K_sub;
            p2.A = Rkron \ p2.A;
            p2.B = Rkron \ p2.B;

            K_sub = p2;
            
            
            
            % K_sub_full = rhotrafo(K_feed_full, 1/rho);
            %connect the internal model: form the controller

            

            K = lft(model, K_sub);
            % K_full = lft(model, K_sub_full);


            % alg_full = lft(obj.sys.P, K_full);


            K_report = struct;            
            K_report.K = K;
            K_report.model = model;
            K_report.K_sub = K_sub;
            K_report.alg_trans = alg_trans;  

        end
               
    end
end