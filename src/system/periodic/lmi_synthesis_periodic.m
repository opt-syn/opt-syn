classdef lmi_synthesis_periodic < lmi_synthesis_interface
    %LMI_SYNTHESIS_PERIODIC synthesisLMIs for algorithmic interconnections
    %involving periodic linear networks and controllers
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
        function obj = lmi_synthesis_periodic(sys,config)
            %LMI_SYNTHESIS_PERIODIC undefined
            %   undefined
            obj@lmi_synthesis_interface(sys, config);
        end

        %% definition of variables and helpers

        function ns = Nss(obj)
            %NSS: Number of subsystems            
            ns = obj.sys.Nss;
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


            


            GX_cell = cell(obj.Nss, 1);
            GY_cell = cell(obj.Nss, 1);
            S_cell = cell(obj.Nss, 1);

            if obj.opts.COMMON
                %common storage function among all subsystems
                
                [GX, GY, cons] = obj.define_storage_G(cons, alg_psi{1}, '');
                n = ssize(GX, 1);
                % vars_diss= struct('GX', GX, 'GY', GY, 'S', eye(n));

                G = vars_diss.G;
                GX_cell = cell(obj.Nss, 1);
                for i = 1:obj.Nss
                    GX_cell{i} = GX;
                    GY_cell{i} = GY;
                    S_cell{i} = eye(n);
                end

            else
                %define a storage function for each subsystem

                for i = 1:obj.Nss
                    [GX, GY, cons] = obj.define_storage_G(cons, alg_psi{i}, num2str(i));
                    n = ssize(GX, 1);
                    GX_cell{i} = GX;
                    GY_cell{i} = GY;
                    S_cell{i} = eye(n);
                end

            end

            vars_diss = struct;
            vars_diss.GX = GX_cell;
            vars_diss.GY = GY_cell;
            vars_diss.S  = S_cell;

        end


        function [vars_K, cons] = create_vars_controller(obj, cons, alg_psi, name)
            %CREATE_VARS_CONTROLLER create the nonlinearly-transformed
            %controller matrices

            %get the dimensions

            vars_K = cell(obj.Nss, 1);

            if nargin < 4
                name = [];
            end

            for i = 1:obj.Nss
                name_curr = [name, '_', num2str(i)];
                alg_curr = alg_psi{i};
                [vars_K{i}, cons] = create_vars_controller@lmi_synthesis_interface(obj, cons, alg_curr, name_curr);
            end
            
        end

        function vars_inv= get_vars_involved(obj, vars, ind)
            %GET_VARS_INVOLVED get variables involved in the current mode

            vars_inv= struct;
            vars_inv.diss.GX = vars.diss.GX{ind};
            vars_inv.diss.GY = vars.diss.GY{ind};
            vars_inv.diss.S  = vars.diss.S{ind};

            vars_inv.reg.Pi = vars.reg.Pi{ind};
            vars_inv.reg.Gam = vars.reg.Gam{ind};
            vars_inv.reg.Phi = vars.reg.Phi{ind};

        end

        function [cons, objective, con_M] = quad(obj, vars, cons, diss)
            %QUAD: certificate of infinite-horizon quadratic performance



            %get the variables of the problem
            vcurr = obj.get_vars_diss(vars, diss.ind_curr);
            vnext = obj.get_vars_diss(vars, diss.ind_curr);
            
            Gcurr = obj.get_storage(vcurr.diss, vcurr.reg);
            Gnext = obj.get_storage(vnext.diss, vnext.reg);




            %system block with {A, B, G}
            sysb = obj.sys_block(diss.plant, Gnext, Gcurr);


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
            cons = obj.con_terminal(Gcurr, cons, diss.iqc_rob);
        end

    end
end