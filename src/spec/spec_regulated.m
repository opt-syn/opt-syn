classdef spec_regulated < spec_interface
    %SPEC_REGULATED enforce  a regulation property 
    %
    %this is used as an abstraction for reduced-order control only
    properties
        type = 'regulated';        
    end

    methods
        function obj = spec_regulated(iwp, izp)
            %SPEC_STABLE Construct an instance of this class
            %   Detailed explanation goes here

            if nargin == 0
                rho = 1;
            end
            obj@spec_interface(iwp, izp);            
        end       

    end
end

