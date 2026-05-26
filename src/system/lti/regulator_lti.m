classdef regulator_lti < regulator_interface
    %REGULATOR_LTI Regulator for LTI systems
    

    
    methods
        function obj = regulator_lti(sys)
            %REGULATOR_LTI Construct an instance of this class
            %   Detailed explanation goes here
            obj@regulator_interface(sys)
        end
        
%% Regulation check and internal models

        %solve the regulator equation for the network (Pi, Gamma, Phi)
        %   used for verification and internal model control
        %   find a parameterization of solutions of the regulator equation

        %construct the internal model

        %solve the regulator equation for the controller (Theta) given
        %(Gamma, Phi): this is verification
        

        function sys = get_model(obj, vars_reg)
            %get_model
            %fetch the internal model (nominal)
            %
            %
            %with edits: allow for selection of model within feasible set

            %TODO: allow for parameterizations based on the variables
            
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
        
    end
end

