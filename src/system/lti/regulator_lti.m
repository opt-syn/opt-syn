classdef regulator_lti < regulator_interface
    %REGULATOR_LTI Regulator for LTI systems
    

    
    methods
        function obj = regulator_lti(sys)
            %REGULATOR_LTI Construct an instance of this class
            %   Detailed explanation goes here
            obj@regulator_interface(sys)
        end
        
        function sys = get_model(obj, vars_reg)
            %get_model
            %fetch the internal model (nominal)
            %
            %
            %with edits: allow for selection of model within feasible set

            %TODO: allow for eterizations based on the variables
            
            if nargin < 2
                Phi = obj.Phi;
                Gam = obj.Gam;
            else
                Phi = vars_reg.Phi;
                Gam = vars_reg.Gam;
            end

            sys =obj.fetch_model(obj.S, Phi, Gam);
            

        end

        function plant_model = connect_model(obj, plant, rho)
            %connect the model (nominal regulator equation)

            if nargin < 3
                rho = 1;
            end

            model = obj.get_model();
            model_rho = rhotrafo(model, rho);
            plant_model = lft(plant, model_rho);
        end

        %% methods for reduced-order control

        function alg_aug = sys_regulated_aug(obj)
            %SYS_REGULATED_AUG augment the system by the regulated
            %disturbance

            % if nargin < 3
            %      =[];
            % end

            [A, B1, B2, C1, D11, D12, C2, D21, D22] = obj.sys.ss_zy_wu();
            [Bd, Ded, Dyd] = obj.d_influence();

            
            Bwp = obj.sys.P.Bwp();
            Czp = obj.sys.P.Czp();


            B_aug = [B1, Bwp, Bd, B2];
            C_aug = [C1; C1; Czp; C2];

            Dzwp = obj.sys.P.Dzwp();
            Dzpw = obj.sys.P.Dzpw();
            Dzpwp = obj.sys.P.Dzpwp();
            Dywp = obj.sys.P.Dywp();
            Dzpu = obj.sys.P.Dzpu();
            nzp = obj.sys.P.nzp;
            nd = size(Bd, 1);
            

            D_aug = [D11, Dzwp, Ded, D12;
                     Dzpw, Dzpwp, zeros(nzp, nd), Dzpu;
                     D11, Dzwp, Ded, D12;
                     D21, Dywp, Dyd, D22];

            alg_aug_P = ss(A, B_aug, C_aug, D_aug, 1);

            nn = obj.sys.P.dump_dim();
            nn.nwp = nn.nwp + size(Bwp, 2);
            nn.nzp = nn.nzp + size(Ded, 1);

            alg_aug = genplant(alg_aug_P, nn);
        end


        
    end
end

