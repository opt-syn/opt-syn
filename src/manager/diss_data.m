classdef diss_data
    %DISS_DATA Container for a dissipation inequality
    %   processed in lmi_analysis and lmi_synthesis routines
    %   formed in the index_specs routine;
    
    properties
        iqc_rob; %(iqc) iqcs from the operator uncertainty
        iqc_data; %(cell) iqcs of the individual operators
        plant;   %(genplant) the plant to analyze/control
        plant_reg; %(genplant) the plant with error channels, used for reduced-ordeer control
        spec;     %(spec) performance specification
        rho=1;    %linear convergence rate
        ndiss = 1;%number of specifications
    end
    
    methods
        function obj = diss_data()
            %construct a dissipation inequality            
        end

    end
end

