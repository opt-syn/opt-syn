classdef lmi_analysis_lti < lmi_analysis_interface
    %LMI_ANALYSIS_LTI analysis LMIs for algorithmic interconnections
    %involving linear-time-invariant (LTI) networks and controllers
    
    % w(k) \in F(z(k))
    %
    % [x(k+1)] = [A    B     Bp   ][x(k)]   state transition
    % [z(k)  ] = [Cz   Dzw   Dzwp ][w(k)]   output to oracle
    % [zp(k) ] = [Czp  Dzpw  Dzpwp][wp(k)]  output to performance
    %
    % performance specification: wp -> zp from (spec)
    %
    %   Implemented
    %       stability
    %       e2e
    %       quad
    %       p2p ?
    %
    %   TODO:
    %       h2      
    %       e2p
    %       
    %
    
    methods
        function obj = lmi_analysis_lti(sys, config)
            %LMI_ANALYSIS_LTI Constructor
            obj@lmi_analysis_interface(sys, config);
        end       
        
        
        %% definition of variables
        function [vars_diss, cons]= create_vars_storage(obj, cons, alg_psi, name)
            %create_vars_storage create variables for the dissipation
            %constraints
            %
            %Args:
            %   cons:       accumulated constraints
            %   alg_psi:    the filtered algorithmic interconnection
            %   name:       a name for the variable
            %Returns:
            %   vars_diss:   variables of the problem in the dissipation constraints
            %   cons: accumulated constraints


            if nargin < 4
                name = [];
            end

            [G, cons] = obj.define_storage_G(cons, alg_psi, name);
            vars_diss= struct('G', G);

        end



        %% performance spceifications

        % quadratic performance  
        function [cons, objective, con_M] = quad(obj, vars, cons, diss)
            %QUAD: certificate of quadratic performance
            %
            %Args:
            %   vars:   variables of the problem        
            %   cons:   accumulated constraints
            %   diss (diss_data):   structure describing the dissipation constraint
            %Returns:
            %   cons:   accumulated constraints
            %   objective:  term to be minimized            
            %   con_M:      PSD blocks for the dynamics constraint

     

            G = vars.diss.G;

            %system block with {A, B, G}
            sysb = obj.sys_block(diss.plant, G, G, diss.rho);

            %supply and quadratic performance constraints
            %index the robustness and performance channels separately
            [plant_rob, plant_perf] = obj.partition_perf(diss);


            %supply block with {C, D, M}
            suppb = -obj.supply_block(plant_rob, diss.iqc_rob.M);

            con_M = -(sysb + suppb);           
            
            %quadratic performance block by Schur complement
            [con_M, objective] = obj.quad_performance_augment(diss, vars, con_M, plant_perf);
                        

            %wrap it all up
            sM = ssize(con_M,1);
            cons = append_lmi(cons, con_M - obj.config.tol.M*eye(sM), obj.config.LMILAB); 

            %impose sign constraint
            cons = obj.con_terminal(G, cons, diss.iqc_rob);
        end    

        function [cons, objective, con_M] = h2(obj, vars, cons, diss)
            %H2: certificate of stochastic performance
            %
            %Args:
            %   vars:   variables of the problem        
            %   cons:   accumulated constraints
            %   diss (diss_data):   structure describing the dissipation constraint
            %Returns:
            %   cons:   accumulated constraints
            %   objective:  term to be minimized            
            %   con_M:      PSD blocks for the dynamics constraint


            %separate the performance inputs from the plant inputs

            [plant_rob, plant_perf] = obj.partition_perf_zp(diss);
            

            %pose quadratic constraint
            G = vars.diss.G;
            %system block with {A, B, G}
            sysb = obj.sys_block(plant_rob, G, G, diss.rho);

            %supply block with {C, D, M} and wp
            M_rob = blkdiag(diss.iqc_rob.M, eye(diss.spec.nwp));
            suppb = -obj.supply_block(plant_rob, M_rob);

            con_M = -(sysb + suppb);           


            %output constraint          
            vars_spec = vars.spec{diss.spec.id};
            Omega = diss.spec.get_cov();
            
            con_Z = h2_block(obj, plant_perf, G, vars_spec, Omega);

            objective = diss.spec.get_objective(vars_spec);



            %wrap it all up
            sM = ssize(con_M,1);            
            cons = append_lmi(cons, con_M - obj.config.tol.M*eye(sM), obj.config.LMILAB); 
            cons = append_lmi(cons, con_Z, obj.config.LMILAB); 

            %impose sign constraint
            cons = obj.con_terminal(G, cons, diss.iqc_rob); 
        end   



        %% Peak-to-Peak norm (at each finite horizon)
        function [cons, objective, con_M] = p2p(obj, vars, cons, diss)
            %p2p: certificate of peak to peak induced norm
            %Args:
            %   vars:   variables of the problem        
            %   cons:   accumulated constraints
            %   diss (diss_data):   structure describing the dissipation constraint
            %Returns:
            %   cons:   accumulated constraints
            %   objective:  term to be minimized            
            %   con_M:      PSD blocks for the dynamics constraint

            

     



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
            
            sysb_1 = obj.sys_block(diss.plant, G, G);
            suppb_1 = obj.supply_block(plant_no_p, M_base);


            con_M_1 = sysb_1 + suppb_1;

            %Block 2: with performance (and terminal constraint)

            sysb_2 = obj.sys_block(diss.plant, X_f, G);


            if diss.spec.target
                %TODO: fix this
                error('Spec P2P with exponential discounting needs some work.')

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

            cons = append_lmi(cons, con_M_1 - obj.config.tol.M*eye(sM1), obj.LMILAB);   
            cons = append_lmi(cons, con_M_2 - obj.config.tol.M*eye(sM2), obj.LMILAB);                         

            con_M = {con_M_1, con_M_2};
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

            %highlight exponential stability

            alg_psi_first = alg_psi{1};

            nwp = alg_psi_first.nwp;
            nzp = alg_psi_first.nzp;

            alg_psi_exp = alg_psi_first(1:(end - (nzp)), 1:(end - (nwp)));

            %evaluate the variables
            alg_psi_exp = rhotrafo(alg_psi_exp, sol.rho);
            iqc_op_all = struct('iqc', diss{1}.iqc_rob);

            %TODO: finish validation
            % sol.gain = obj.validate_recovery_gain(alg_psi_exp, iqc_op_all);

            % sol.G = obj.get_storage(sol.vars.diss, sol.vars.reg);
        end


    end

    
end


