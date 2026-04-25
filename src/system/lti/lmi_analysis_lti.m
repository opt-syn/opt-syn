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
        
        
            %CON_TERMINAL terminal constraint
        



        function [cons, objective, con_M] = quad(obj, vars, cons, diss)
            %QUAD: certificate of quadratic performance

            G = vars.diss.G;

            sysb = obj.sys_block(diss.plant, G, G, diss.spec.rho);


            M_quad = obj.merge_spec_M(diss.iqc_rob, diss.spec);
            suppb = obj.supply_block(diss.plant, M_quad);


            %wrap it all together
            objective = 0;

            con_M = sysb + suppb;


            sM = ssize(con_M,1);
            cons = append_lmi(cons, con_M - obj.tol.M*eye(sM), obj.LMILAB); 

        end





        function [cons, objective, con_M] = e2e_target(obj, vars, cons, diss)
            %E2E_TARGET: use a Schur complement to minimize the energy to
            %energy gain of the transfer function

            G = vars.diss.G;

            
            sysb = obj.sys_block(diss.plant, G, G, diss.spec.rho);

            %variable to optimize
            mu = vars.spec{diss.spec.id}.mu_l2;

            %separate the performance outputs   
            nzp = length(diss.spec.izp);
            ind_sep = (diss.iqc_rob.np) + (1:nzp);
            nz = ssize(diss.plant.D, 1);

            ind_diff_sep = setdiff(1:nz, ind_sep);
            Iz = eye(nz);
            Izp = Iz(ind_diff_sep, :);

            %the plant without the schur-complemented-out performance input
            alg_screen = Izp*diss.plant;


            Ezp = sparse(1:length(ind_sep), ind_sep, ones(length(ind_sep)), ...
                length(ind_sep), ssize(diss.plant.C, 1));


            CDp = Ezp * [diss.plant.C, diss.plant.D];
                       

            %form the supply
            nwp = length(diss.spec.iwp);
            M_base = blkdiag(diss.iqc_rob.M, -mu * eye(nwp));
            
            objective = mu;            
            suppb = obj.supply_block(alg_screen, M_base);


            %wrap it all together           
            con_M_top = sysb + suppb;
            con_M = [con_M_top, CDp'; CDp, mu*eye(nzp)];


            sM = ssize(con_M,1);
            cons = append_lmi(cons, con_M - obj.tol.M*eye(sM), obj.LMILAB);             
        end




    end

    
end

