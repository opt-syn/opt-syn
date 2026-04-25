classdef lmi_analysis_lti < lmi_analysis_interface
    %LMI_ANALYSIS_LTI analysis LMIs for algorithmic interconnections
    %involving linear-time-invariant (LTI) networks and controllers
    %
    % w(k) \in F(z(k))
    %
    % [x(k+1)] = [A    B     Bp   ][x(k)]   state transition
    % [z(k)  ] = [Cz   Dzw   Dzwp ][w(k)]   output to oracle
    % [zp(k) ] = [Czp  Dzpw  Dzpwp][wp(k)]  output to performance
    %
    % performance specification: wp -> zp from (spec)
    %
    %   Implemented
    %
    %
    %   TODO:
    %       stability
    %       h2
    %       e2p
    %       p2p
    %
    
    methods
        function obj = lmi_analysis_lti(sys)
            %LMI_DISPATCH_LTI Construct an instance of this class
            %   Detailed explanation goes here
            obj@lmi_analysis_interface(sys);
        end

        %function [vars_diss, cons]= create_vars_storage(obj, alg_psi, cons, name)
        %is the default
        
        


        %% Quadratic performance (infinite horizon)
        function [cons, objective, con_M] = quad(obj, vars, cons, diss)
            %QUAD: certificate of infinite-horizon quadratic performance

            G = vars.diss.G;

            sysb = obj.sys_block(diss.plant, G, G, diss.spec.rho);


            M_quad = obj.merge_spec_M(diss.iqc_rob, diss.spec);
            suppb = obj.supply_block(diss.plant, M_quad);


            %wrap it all together
            objective = 0;

            con_M = sysb + suppb;


            sM = ssize(con_M,1);
            cons = append_lmi(cons, con_M - obj.tol.M*eye(sM), obj.LMILAB); 

            %impose sign constraint
            cons = obj.con_terminal(G, cons, diss.iqc_rob);
        end        



        function [cons, objective, con_M] = e2e_target(obj, vars, cons, diss)
            %E2E_TARGET: use a Schur complement to minimize the energy to
            %energy gain of the transfer function

            G = vars.diss.G;

            
           
            sysb = obj.sys_block(diss.plant, G, G, diss.spec.rho);

            %variable to optimize
            mu = vars.spec{diss.spec.id}.mu_l2;

            [plant_no_p, CDp] = obj.separate_performance_output(diss);

            %form the supply
            nwp = length(diss.spec.iwp);
            M_base = blkdiag(diss.iqc_rob.M, -mu * eye(nwp));
            
            objective = mu;            
            suppb = obj.supply_block(plant_no_p, M_base);


            %wrap it all together           
            con_M_corner = sysb + suppb;
            nzp = ssize(CDp, 1);
            con_M = [con_M_corner, CDp'; CDp, mu*eye(nzp)];


            sM = ssize(con_M,1);
            cons = append_lmi(cons, con_M - obj.tol.M*eye(sM), obj.LMILAB);   
            
            
            %impose sign constraint
            cons = obj.con_terminal(G, cons, diss.iqc_rob);
        end

        %% Peak-to-Peak norm (at each finite horizon)

        function [cons, objective, con_M] = p2p(obj, vars, cons, diss)
            %p2p: certificate of peak to peak induced norm
            %
            % sup norm(zp, 2) / norm(wp, 2) <= objective

            % verification by Theorem 4 of https://www.sciencedirect.com/science/article/pii/S2405896323008194

            %storage matrix
            G = vars.diss.G;
                      
            %terminal constraint
            X = diss.iqc_rob.X;
            nf = ssize(X);
            n = ssize(G, 1);
            Ef = [eye(nf); zeros(n-nf, nf)];

            X_f = Ef * X * Ef';

            

            %variables in the problem
            vars_spec = vars.spec{diss.spec.id};
            mu = vars_spec.mu_p2p;            

            %form the plant
            [plant_no_p, CDp] = obj.separate_performance_output(diss);


            %Block 1: without performance
            
            nwp = length(diss.spec.iwp);
            M_base = blkdiag(diss.iqc_rob.M, -mu * eye(nwp));
            
            sysb_1 = obj.sys_block(diss.plant, G, G, diss.spec.rho);
            suppb_1 = obj.supply_block(plant_no_p, M_base);


            con_M_1 = sysb_1 + suppb_1;

            %Block 2: with performance (and terminal constraint)

            sysb_2 = obj.sys_block(diss.plant, X_f, G, diss.spec.rho);


            if diss.spec.target
                

                %optimize over the gain

                %hardcode the supply function

                rho = diss.spec.rho;
                rrecip = rho/(1-rho);
                gam = vars_spec.gam_p2p;
                M_u = -eye(diss.spec.nwp) * rrecip * (gam - mu);                

                

                
                M_p2p = blkdiag(diss.iqc_rob.M, M_u);
                suppb_2 = obj.supply_block(plant_no_p, M_p2p);

                nzp = ssize(CDp, 1);

                %reciprocal by Schur complement
                M_yr =  (gam)*eye(nzp) *  (rrecip)^(-1);                

                con_M_2_corner = suppb_2 + sysb_2;

                %expand out the Schur complement
                con_M_2 = [con_M_2_corner, CDp'; CDp, M_yr];

                objective = gam;
            else
                
            
                
    
                %remember to take proper block diagonals and indexes!
                M_p2p = obj.merge_spec_M(diss.iqc_rob, diss.spec, vars_spec);
                
                
                suppb_2 = obj.supply_block(diss.plant, M_p2p);
                
                con_M_2 = sysb_2 + suppb_2;

                objective = 0;
            end

            %wrap it all together           

            sM1 = ssize(con_M_1,1);  sM2 = ssize(con_M_2,1);

            cons = append_lmi(cons, con_M_1 - obj.tol.M*eye(sM1), obj.LMILAB);   
            cons = append_lmi(cons, con_M_2 - obj.tol.M*eye(sM2), obj.LMILAB);                         

            con_M = {con_M_1, con_M_2};
        end

    end

    
end


