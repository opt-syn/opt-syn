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

        function [vars, cons, objective, con_M] = cons_dynamic(obj, vars, cons, diss)
            %CONS_DYNAMIC form the dissipation and sign constraints


            %go to the rotating coordinate frame
            diss.plant = obj.rotate_plant(diss.plant);

            %now call the cons_dynamic routine for LTI systems
            [vars, cons, objective, con_M] = cons_dynamic@lmi_synthesis_lti(obj, vars, cons, diss);
        end

        function sol = process_recovery(obj, sol, lmi_out, alg_psi, diss)
            %recover the controller

            if nargin < 5
                diss = [];
            end

            %get the system with the internal model
            dissend = struct;
            dissend.plant = obj.rotate_plant(alg_psi, 1);
            dissend.rho = sol.rho;            
            P_trans =  obj.connect_model(dissend);

            %evaluate the variables
            [sol] = obj.recover_subcontroller(alg_psi, P_trans, sol);

            % sol.G = obj.get_storage(sol.vars.diss, sol.vars.reg);
        end


        function plant_rot = rotate_plant(obj, plant, direction)
            
            %rotate_plant: apply the periodic-orbit rotation to the
            %time-varying system, producing an LTI system
            if nargin < 3
                direction = 1;
            end

            c = size(obj.sys.R, 1);
            n = size(plant.A, 1);
            Rkron = kron(eye(n/c), obj.sys.R)^(direction);
            
            plant_rot = plant;
            if isa(plant, 'genplant')
                plant_rot.P.A = Rkron * plant.P.A;
                plant_rot.P.B = Rkron * plant.P.B;
            else
                plant_rot.A = Rkron * plant.A;
                plant_rot.B = Rkron * plant.B;
            end


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


            [sol] = recover_subcontroller@lmi_synthesis_lti(obj, alg_psi, P_trans, sol);


            %revert the coordinate transformation
            
            sol.cert.alg = obj.rotate_plant(sol.cert.alg, -1);            
            sol.cert.model = obj.rotate_plant(sol.cert.model, -1);            
            sol.K= obj.rotate_plant(sol.K, -1);            
            sol.cert.K_sub = obj.rotate_plant(sol.cert.K_sub, -1);                        

        end
               
    end
end
