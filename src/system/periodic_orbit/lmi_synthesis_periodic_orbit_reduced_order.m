classdef lmi_synthesis_periodic_orbit_reduced_order < lmi_synthesis_lti_reduced_order
    %LMI_SYNTHESIS_PERIODIC_ORBIT_REDUCED_ORDER reduced-order synthesis LMIs for algorithmic interconnections
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
        function obj = lmi_synthesis_periodic_orbit_reduced_order(sys,config)
            %LMI_SYNTHESIS_PERIODIC_ORBIT_REDUCED_ORDER constructor
            obj@lmi_synthesis_lti_reduced_order(sys, config);
        end

        function [vars, cons, objective, con_M] = cons_dynamic(obj, vars, cons, diss)
            %CONS_DYNAMIC form the dissipation and sign constraints
            %
            %Args:
            %   vars:   variables of the problem        
            %   cons:   accumulated constraints
            %   diss (diss_data):   structure describing the dissipation constraint
            %
            %Returns:
            %   cons:   accumulated constraints
            %   objective:  term to be minimized            
                      


            %go to the rotating coordinate frame
            % diss.plant = obj.rotate_plant(diss.plant);

            %now call the cons_dynamic routine for LTI systems
            [vars, cons, objective, con_M] = cons_dynamic@lmi_synthesis_lti(obj, vars, cons, diss);
        end

        function sol = process_recovery(obj, sol, lmi_out, alg_psi, diss)
            %recover the controller
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
            % dissend.plant = obj.rotate_plant(alg_psi, 1);
            dissend.plant_reg = diss{1}.plant_reg;
            dissend.rho = sol.rho;      
            dissend.iqc_data = diss{1}.iqc_data;
            P_trans =  obj.connect_model(dissend);

            %evaluate the variables
            [sol] = obj.recover_subcontroller(alg_psi, P_trans, sol);


            % sol.G = obj.get_storage(sol.vars.diss, sol.vars.reg);
        end


        function plant_rot = rotate_plant(obj, plant, direction)            
            %rotate_plant,  apply the periodic-orbit rotation to the
            %time-varying system, producing an LTI system
            %
            %Args:
            %   plant: original system
            %   direction (bool): forwards (true) or backwards (false)
            %Returns:
            %   plant_rot: rotated LTI plant
            if nargin < 3
                direction = 1;
            end

            c = size(obj.sys.M, 1);
            n = size(plant.A, 1);
            Rkron = kron(eye(n/c), obj.sys.M)^(direction);
            
            plant_rot = plant;
            if isa(plant, 'genplant')
                plant_rot.P.A = Rkron * plant.P.A;
                plant_rot.P.B = Rkron * plant.P.B;
            else
                plant_rot.A = Rkron * plant.A;
                plant_rot.B = Rkron * plant.B;
            end


        end

        function [sol] = recover_subcontroller(obj, alg_psi, P_aug, sol)
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
            
            vars_rec = sol.vars;
            rho = sol.rho;

            %get the subcontroller
            [K_nofeed, Gcl, Ycl] = recover_subcontroller_warp(obj, P_aug, vars_rec);


            %get the controller
            model = obj.reg.get_model(vars_rec.reg); 
            % model_orig = obj.rotate_plant(model, -1);    
            
            modelrho = rhotrafo(model, sol.rho); 
            % P_trans = lft(alg_psi, modelrho);
            P_trans = lft(alg_psi, model);
            % K_rot = obj.rotate_plant(K_nofeed, -1);

            K_report = obj.K_alg_report(P_trans, K_nofeed, model);
            
            
            

            %export and validate
            sol.cert.K=  obj.rotate_plant(K_report.K, -1);
            sol.cert.alg_trans =  K_report.alg_trans;
            sol.cert.alg = lft(obj.sys.P, sol.cert.K);
            sol.cert.model = obj.rotate_plant(model, -1);           
            
            sol.cert.K_sub = obj.rotate_plant(K_report.K_sub, -1);
            alg_rho = rhotrafo(sol.cert.alg_trans, sol.rho);
            sol.gain = obj.validate_recovery_gain(alg_rho, sol.cert.iqc_op_all);

            sol.cert.Gcl = Gcl;
            sol.cert.Ycl = Ycl;
 
        end
               
    end
end
