classdef diss_data
    %DISS_DATA Container for a dissipation inequality
    %   processed in lmi_analysis and lmi_synthesis routines
    
    properties
        iqc_rob; %(cell) iqcs from the uncertainty (robust)
    end
    
    methods
        function obj = diss_data;
            %DISS construct a dissipation inequality
            %   Detailed explanation goes here            
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

