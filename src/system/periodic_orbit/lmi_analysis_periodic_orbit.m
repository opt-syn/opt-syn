classdef lmi_analysis_periodic_orbit < lmi_analysis_periodic
    %LMI_ANALYSIS_PERIODIC_ORBIT analysis LMIs for algorithmic interconnections
    %involving periodic linear networks and controllers
    %
    % w(k) \in F(z(k))
    %
    % [x(k+1)] = [A(k)    B(k)     Bp(k)   ][x(k)]   state transition
    % [z(k)  ] = [Cz(k)   Dzw(k)   Dzwp(k) ][w(k)]   output to oracle
    % [zp(k) ] = [Czp(k)  Dzpw(k)  Dzpwp(k)][wp(k)]  output to performance
    %
    %A(k) = A(k+T) for some known time T
    %furthermore, matrices [R, W] are known with 
    %A(k) = (R^k)' A(0) R^k,  R^T = I (and the same for other channels.
    %
    % this is a specialization of general periodic algorithms
    %
    %
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
        R; %periodicity in the state in dynamics
    end

    methods
        function obj = lmi_analysis_periodic_orbit(sys, config)
            %LMI_DISPATCH_LTI Construct an instance of this class
            %   Detailed explanation goes here
            obj@lmi_analysis_periodic(sys, config);
            obj.R = sys.R;
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
    
            diss.plant = diss.plant{1};
          [cons, objective, con_M] = cons_dynamic@lmi_dispatch_interface(obj, vars, cons, diss);
        end


        %% LMI definitions
        function [cons, objective, con_M] = quad(obj, vars, cons, diss)
            %QUAD: certificate of infinite-horizon quadratic performance


            
            G = vars.diss.G;

            nf = ssize(diss.plant.A,1);

            Rkron = kron(eye(nf/size(obj.R, 1)), obj.R);

            Gnew = Rkron'*G*Rkron;

            %system block with {A, B, G}
            sysb = obj.sys_block(diss.plant, Gnew, G);


            %supply block with {C, D, M}
            vars_spec = vars.spec{diss.spec.id};
            M_quad = obj.merge_spec_M(diss.iqc_rob, diss.spec, vars_spec);
            suppb = obj.supply_block(diss.plant, M_quad);


            %wrap it all together
            objective = 0;

            con_M = sysb + suppb;


            sM = ssize(con_M,1);
            cons = append_lmi(cons, con_M - obj.config.tol.M*eye(sM), obj.config.LMILAB); 

            %impose sign constraint
            cons = obj.con_terminal(G, cons, diss.iqc_rob);

        end        

        function [cons, objective, con_M] = e2e_target(obj, vars, cons, diss)
            %E2E_TARGET: use a Schur complement to minimize the energy to
            %energy gain of the transfer function

            G = vars.diss.G;
                        
            nf = ssize(diss.plant.A,1);

            Rkron = kron(eye(nf/size(obj.R, 1)), obj.R);

            Gnew = Rkron'*G*Rkron;
           
            sysb = obj.sys_block(diss.plant, G, G);

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
            cons = obj.con_terminal(G, cons, diss.iqc_rob);

        end



        %% helper routines
        function ns = Nss(obj)
            %NSS: Number of subsystems            
            ns = obj.sys.Nss;
        end



        function [vars_diss, cons]= create_vars_storage(obj, cons, alg_psi, name)
            %create_vars_storage create variables for the dissipation
            %constraints
            %
            %Input:
            %   cons:       accumulated constraints
            %   alg_psi:    the filtered algorithmic interconnection
            %   name:       a name for the variable

            if nargin < 4
                name = [];
            end

            [G, cons] = obj.define_storage_G(cons, alg_psi, name);
            vars_diss= struct('G', G);

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