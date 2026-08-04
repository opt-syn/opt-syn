classdef lmi_synthesis_periodic_orbit < lmi_synthesis_lti
    %LMI_SYNTHESIS_PERIODIC_ORBIT synthesis LMIs for algorithmic interconnections
    %involving periodic linear networks and controllers
    
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
            %LMI_SYNTHESIS_PERIODIC constructor
            %   undefined
            obj@lmi_synthesis_lti(sys, config);
        end
        function sol = process_recovery(obj, sol, lmi_out, alg_psi, diss)
            %recover the controller
            %
            %Args:
            %   sol: solution structure
            %   lmi_out: output from solver
            %   alg_psi:   the filtered algorithmic interconnection
            %   diss (diss_data):   structure describing the dissipation constraint            
            %
            %Returns:  
            %   sol: solution structure

            if nargin < 5
                diss = [];
            end

            %get the system with the internal model
            dissend = struct;
            dissend.plant = alg_psi;
            dissend.rho = sol.rho;            
            P_trans =  obj.connect_model(dissend);

            %evaluate the variables
            [sol] = obj.recover_subcontroller(alg_psi, P_trans, sol);

            % sol.G = obj.get_storage(sol.vars.diss, sol.vars.reg);
        end

        function [sol] = recover_subcontroller(obj, alg_psi, P_trans, sol)
            %RECOVER_SUBCONTROLLER recover the subcontroller of the current
            %mode/control
            %
            %Args:
            %   alg_psi:   the filtered algorithmic interconnection
            %   P_trans:    the transformed generalized plant before IQC
            %   sol: solution structure
            %
            %Returns:
            %   sol: solution structure
            

            [sol] = recover_subcontroller@lmi_synthesis_lti(obj, alg_psi, P_trans, sol);


            %revert the coordinate transformation
            
            sol.cert.alg = obj.sys.rotate_plant(sol.cert.alg, -1);            
            sol.cert.model = obj.sys.rotate_plant(sol.cert.model, -1);            
            sol.cert.K= obj.sys.rotate_plant(sol.cert.K, -1);            
            sol.cert.K_sub = obj.sys.rotate_plant(sol.cert.K_sub, -1);                        

        end
               
    end
end
