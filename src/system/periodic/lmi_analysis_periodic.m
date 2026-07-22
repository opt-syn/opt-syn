classdef lmi_analysis_periodic < lmi_analysis_interface
    %LMI_ANALYSIS_PERIODIC analysis LMIs for algorithmic interconnections
    %involving periodic linear networks and controllers
    
    % w(k) \in F(z(k))
    %
    % [x(k+1)] = [A(k)    Bw(k)     Bwp(k)  ][x(k)]   state transition
    % [z(k)  ] = [Cz(k)   Dzw(k)   Dzwp(k) ][w(k)]   output to oracle
    % [zp(k) ] = [Czp(k)  Dzpw(k)  Dzpwp(k)][wp(k)]  output to performance
    %
    %A(k) = A(k+T) for some known time T
    %
    %instances of these algorithms include cyclic coordinate descent
    %methods. Periodic systems can also be unrolled into an LTI system
    %(monodromy methods): a single large LMI system rather than multiple 
    % coupled smaller LMI systems
    %
    %   Implemented
    %       stability
    %       e2e
    %
    %   TODO:
    %       h2      
    %       e2p
    %       
    %


    properties
        opts = struct("COMMON", false);
    end

    methods
        function obj = lmi_analysis_periodic(sys, config)
            %LMI_DISPATCH_LTI Constructor
            obj@lmi_analysis_interface(sys, config);
        end


        %% main call
        function [vars, cons, objective] = cons_dynamic(obj, vars, cons, diss)
            %CONS form the dissipation and sign constraints
            %
            %Args:
            %   vars:   variables of the problem        
            %   cons:   accumulated constraints
            %   diss (diss_data):   structure describing the dissipation
            %       constraint
            %
            %Returns:
            %   cons:   accumulated constraints
            %   objective:  term to be minimized            
                      
            
            %need to look up the right constraint            

            %Upper-levels: iterate over the systems
            objective = 0;
            for i = 1:obj.Nss
                %extract the information of subsystem i
                diss_curr = diss;
                diss_curr.plant = diss.plant{i};
                diss_curr.ind_curr = i;
                diss_curr.ind_next = 1+mod(i, obj.Nss);
                

                [cons, objective_curr, con_M] = obj.con_dynamic_single(vars, cons, diss_curr);
                
                
                %TODO: take the max over the different subsystems
                %but the same objective is sent to each subsystem, so it's
                %all the same? Check this
                if i==1
                    objective = objective + objective_curr;
                end
            end         
                      
        end


        %% LMI definitions
        function [cons, objective, con_M] = quad(obj, vars, cons, diss)
            %QUAD: certificate of infinite-horizon quadratic performance
            %
            %Args:
            %   vars:   variables of the problem        
            %   cons:   accumulated constraints
            %   diss (diss_data):   structure describing the dissipation
            %       constraint
            %
            %Returns:
            %   cons:   accumulated constraints
            %   objective:  term to be minimized            
                      

            
            
            Gcurr = vars.diss.G{diss.ind_curr};
            Gnext = vars.diss.G{diss.ind_next};
            
            
            %allow for differing one-step exponential growths along arcs
            %graph Lyapunov function format
            
            
            %system block with {A, B, G}
            sysb = obj.sys_block(diss.plant, Gnext, Gcurr);


            %supply block with {C, D, M}
            vars_spec = vars.spec{diss.spec.id};
            M_quad = obj.merge_spec_M(diss.iqc_rob, diss.spec, vars_spec);
            suppb = obj.supply_block(diss.plant, M_quad);


            %wrap it all together
            objective = 0;

            con_M = -(sysb + suppb);


            sM = ssize(con_M,1);
            cons = append_lmi(cons, con_M - obj.config.tol.M*eye(sM), obj.config.LMILAB); 

            %impose sign constraint
            cons = obj.con_terminal(Gcurr, cons, diss.iqc_rob);
        end        

        function [cons, objective, con_M] = e2e_target(obj, vars, cons, diss)
            %E2E_TARGET: use a Schur complement to minimize the energy to
            %energy gain of the transfer function
            %Args:                   
            %   cons:   accumulated constraints
            %   specs: performance specifications
            %
            %Returns:            
            %  vars_spec:   variables for performance specification
            %   cons:   accumulated constraints
            
            Gcurr = vars.diss.G{diss.ind_curr};
            Gnext = vars.diss.G{diss.ind_next};



            
           
            sysb = obj.sys_block(diss.plant, Gnext, Gcurr);

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
            cons = append_lmi(cons, con_M - obj.config.tol.M*eye(sM), obj.config.LMILAB);   
            
            
            %impose sign constraint
            cons = obj.con_terminal(Gcurr, cons, diss.iqc_rob);
        end



        %% helper routines
        function ns = Nss(obj)
            %NSS: Number of subsystems            
            ns = obj.sys.Nss;
        end

        function [vars_diss, cons]= create_vars_storage(obj, cons, alg_psi, name)
            %create_vars_storage create variables for the dissipation
            %constraints. One for each subsystem
            %
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

            G_cell = cell(obj.Nss, 1);

            if obj.config.switched.common
                %common storage function among all subsystems
                [vars_diss, cons] = create_vars_storage@lmi_analysis_interface(obj, alg_psi, cons, name);

                G = vars_diss.G;
                G_cell = cell(obj.Nss, 1);
                for i = 1:obj.Nss
                    G_cell{i} = G;
                end
                
            else
                %define a storage function for each subsystem
                for i = 1:obj.Nss
                    [G_curr, cons] = obj.define_storage_G(cons, alg_psi{i}, num2str(i));
                    G_cell{i} = G_curr;
                end

            end

            vars_diss = struct;
            vars_diss.G = G_cell;

        end

        function [vars_spec, cons] = create_vars_spec(obj, cons, specs)
            %CREATE_VARS_SPEC declare variables for the specifications
            %Args:                   
            %   cons:   accumulated constraints
            %   specs: performance specifications
            %
            %Returns:            
            %  vars_spec:   variables for performance specification
            %   cons:   accumulated constraints
            
            
            %right now the variables are in the (spec) object.
            nspec = length(specs);

            vars_spec = cell(nspec, 1);            
            for i = 1:nspec                
                    name = [];
                    [vars_spec{i}, cons] = specs{i}.create_vars(cons, name, obj.config);                
            end
        end
    end
end