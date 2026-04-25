classdef lmi_analysis_interface < lmi_dispatch_interface
    %LMI_ANALYSIS_INTERFACE 
    %Linear Matrix Inequality constraints for analysis of algorithmic
    %interconnections.
    %
    %
    %this is overridden by specialized analysis routines for system types:
    %   lti
    %   periodic
    %   switched robust
    %   switched jump
    
    properties
        Property1
    end
    
    methods
        function obj = lmi_analysis_interface(sys)
            %LMI_ANALYSIS_INTERFACE Construct an instance of this class
            %   Detailed explanation goes here
            obj@lmi_dispatch_interface(sys);
        end

    end
end

