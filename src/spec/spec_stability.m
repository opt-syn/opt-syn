classdef spec_stability < spec_interface
    %SPEC_STABLE enforce exponential stability of the algorithmic
    %interconnection
    
    properties
        type = 'stability';
        rho = 1;
    end
    
    methods
        function obj = spec_stability(rho)
            %SPEC_STABLE Construct an instance of this class
            %   Detailed explanation goes here

            if nargin == 0
                rho = 1;
            end
            obj@spec_interface([], []);

            obj.rho = rho;
        end       

    end
end

