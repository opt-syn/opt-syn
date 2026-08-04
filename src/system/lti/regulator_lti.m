classdef regulator_lti < regulator_interface
    %REGULATOR_LTI Regulator for LTI systems    
    
    methods
        function obj = regulator_lti(sys)
            %REGULATOR_LTI Constructor 
            obj@regulator_interface(sys)
        end
        
        function model = get_model(obj, vars_reg)
            %GET_MODEL fetch the internal model (nominal)
            %
            %Args:
            %   vars_reg:   variables of the problem        (regulator)            
            %Return:
            %   model: the full-order internal model
            
            
            
            if nargin < 2
                Phi = obj.Phi;
                Gam = obj.Gam;
            else
                Phi = vars_reg.Phi;
                Gam = vars_reg.Gam;
            end

            model =obj.fetch_model(obj.S, Phi, Gam);
            

        end

        function plant_model = connect_model(obj, plant, rho)
            %connect the plant to the model (nominal regulator equation)
            %and discount by rho
            %
            %Args:
            %   plant: original system
            %   rho:    exponential weighting
            %Return:
            %   plant_model: plant and model together

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
            %disturbance. ONLY used for reduced-order control
            %
            %Return:
            %   alg_aug: augmented plant with output channels [z, zp, e, y] and input channels [w, wp, d, u]
            %

            

            [A, B1, B2, C1, D11, D12, C2, D21, D22] = obj.ss_zy_wu();
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
                     Dzpw, Dzpwp, zeros(nzp, size(Ded, 2)), Dzpu;
                     D11, Dzwp, Ded, D12;
                     D21, Dywp, Dyd, D22];

            alg_aug_P = ss(A, B_aug, C_aug, D_aug, 1);

            nn = obj.sys.P.dump_dim();


            %place together w and wp, separate d into a new channel group
            %(wp)
            nn.nw = nn.nw + nn.nwp;
            nn.nz = nn.nz + nn.nzp;
            nn.nwp = size(Bd, 2);
            nn.nzp = size(Ded, 1);

            alg_aug = genplant(alg_aug_P, nn);
        end


        
    end
end

