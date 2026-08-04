classdef plant_dims
    %PLANT_DIMS track the dimensions used in the generalized plant.
    %Partition the state space system accordingly.
    
    properties
        s = 0; %number of operators
        nw = 0; %number of operator outputs
        nwp = 0; %number of performance inputs
        nu = 0; %number of controller outputs
        nz = 0; %number of operator inputs
        nzp = 0; %number of performance outputs
        ny = 0; %number of controller inputs

    end
    
    methods
        function obj = plant_dims(nw, nwp, nu, nz, nzp, ny, s)
            %PLANT_DIMS Constructor

        end

    end
end

