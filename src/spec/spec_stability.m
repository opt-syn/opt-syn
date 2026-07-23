classdef spec_stability < spec_interface
    %SPEC_STABLE enforce linear convergence/exponential stability 
    % of the algorithmic interconnection     
    
    methods
        function obj = spec_stability(rho)
            %SPEC_STABLE Constructor
            if nargin == 0
                rho = 1;
            end
            obj@spec_interface([], []);
            obj.type = 'stability';     
            obj.rho = rho;
        end       

    end
end

