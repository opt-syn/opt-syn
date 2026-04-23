classdef spec_stability < spec_interface
    %SPEC_STABLE enforce exponential stability of the algorithmic
    %interconnection
    
    properties
        type = 'stab';
    end
    
    methods
        function obj = spec_stability(rho)
            %SPEC_STABLE Construct an instance of this class
            %   Detailed explanation goes here

            if nargin > 2
                rho = 1;
            end

            obj@spec_interface([], [], rho);
        end       

    end
end

