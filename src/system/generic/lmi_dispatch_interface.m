classdef lmi_dispatch_interface
    %LMI_DISPATCH_INTERFACE analysis and synthesis LMIs for the
    %algorithmic interconnections
    %
    %This contains generic routines common among both analysis and
    %synthesis for every system type
    %  
    properties
        sys;
        LMILAB = 1;
        tol = struct('M', 1e-7, ... %tolerance for dissipation constraints
            'X', 1e-7, ...          %tolerance for sign/terminal cost constraints 
            'G_max', 100)           %upper bound on norm of storage matrix
    end
    
    methods
        function obj = lmi_dispatch_interface(sys)
            %LMI_DISPATCH_INTERFACE Construct an instance of this class
            %   Detailed explanation goes here
            obj.sys = sys;
        end
        
        
        %KYP lemma terms, commonly found matrices
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

