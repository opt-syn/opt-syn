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
        end




    end
end

