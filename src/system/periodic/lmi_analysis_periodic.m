classdef lmi_analysis_periodic < lmi_analysis_interface
    %LMI_ANALYSIS_PERIODIC analysis LMIs for algorithmic interconnections
    %involving periodic linear networks and controllers
    %
    % w(k) \in F(z(k))
    %
    % [x(k+1)] = [A(k)    B(k)     Bp(k)   ][x(k)]   state transition
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
    %       quad
    %       p2p
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
        function obj = lmi_analysis_periodic(sys)
            %LMI_DISPATCH_LTI Construct an instance of this class
            %   Detailed explanation goes here
            obj@lmi_analysis_interface(sys);
        end


        %% main call
        function [cons, objective, con_M] = cons_dynamic(obj, vars, cons, diss)
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
            %       ind_curr:  the index of the current subsystem
            %       ind_next:  the index of the next subsystem
            %
            %Output:
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


            
            
            Gcurr = vars.diss.G{diss.ind_curr};
            Gnext = vars.diss.G{diss.ind_next};


            
            
            %allow for differing one-step exponential growths along arcs
            %graph Lyapunov function format
            rho = obj.get_rho(diss.spec.rho, diss.ind_curr);
            
            %system block with {A, B, G}
            sysb = obj.sys_block(diss.plant, Gnext, Gcurr, rho);


            %supply block with {C, D, M}
            vars_spec = vars.spec{diss.spec.id};
            M_quad = obj.merge_spec_M(diss.iqc_rob, diss.spec, vars_spec);
            suppb = obj.supply_block(diss.plant, M_quad);


            %wrap it all together
            objective = 0;

            con_M = sysb + suppb;


            sM = ssize(con_M,1);
            cons = append_lmi(cons, con_M - obj.tol.M*eye(sM), obj.LMILAB); 

            %impose sign constraint
            cons = obj.con_terminal(Gcurr, cons, diss.iqc_rob);
        end        

        function [cons, objective, con_M] = e2e_target(obj, vars, cons, diss)
            %E2E_TARGET: use a Schur complement to minimize the energy to
            %energy gain of the transfer function

            Gcurr = vars.diss.G{diss.ind_curr};
            Gnext = vars.diss.G{diss.ind_next};



            rho = obj.get_rho(diss.spec.rho, diss.ind_curr); 
           
            sysb = obj.sys_block(diss.plant, Gnext, Gcurr, rho);

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
            cons = obj.con_terminal(Gcurr, cons, diss.iqc_rob);
        end



        %% helper routines
        function ns = Nss(obj)
            %NSS: Number of subsystems            
            ns = obj.sys.Nss;
        end

        function rho = get_rho(obj, rho_list, ind_arc)
            %GET_RHO: get the exponential growth parameter
            if length(rho_list) > 1
                rho = rho_list(ind_arc);
            else
                rho = rho_list;
            end
        end

        function [vars_diss, cons]= create_vars_storage(obj, cons, alg_psi, name)
            %create_vars_storage create variables for the dissipation
            %constraints. One for each subsystem
            %
            %
            %a cell of G(s) functions

            if nargin < 4
                name = [];
            end

            G_cell = cell(obj.Nss, 1);

            if obj.opts.COMMON
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
                    G_curr = obj.define_storage_G(cons, alg_psi{i}, num2str(i));
                    G_cell{i} = G_curr;
                end

            end

            vars_diss = struct;
            vars_diss.G = G_cell;

        end

        function [vars_spec, cons] = create_vars_spec(obj, cons, specs)
            %CREATE_VARS_SPEC declare variables for the specifications

            %maybe put this somewhere else?
            %
            %right now the variables are in the (spec) object.
            nspec = length(specs);

            vars_spec = cell(nspec, 1);
            % vars_spec = cell(nspec, obj.Nss);
            for i = 1:nspec
                % for j = 1:obj.Nss
                    % name = ['_', num2str(j)];
                    name = [];
                    [vars_spec{i}, cons] = specs{i}.create_vars(cons, name);
                % end
            end
        end
    end
end