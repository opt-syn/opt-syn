classdef reg_cl_out
    %REG_CL_OUT Output structure for closed-loop regulator equation
    
    properties
        S; %exosystem signal generator
        R; %exosystem output to plant
        Pi; %regulator equation solution, tracking state of network 
        Gam; %regulator equation solution, tracking input of controller 
        Phi; %regulator equation solution, tracking output of controller 
        Th; %regulator equation solution, tracking state of controller

    end
    
    methods
        function obj = reg_cl_out()
            %REG_CL_OUT Constructor            
        end
        

    end
end

