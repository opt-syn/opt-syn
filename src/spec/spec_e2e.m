classdef spec_e2e < spec_interface
    %SPEC_E2E Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        type='e2e';
        gain = 1;
    end
    
    methods
        function obj = spec_e2e(GAIN, iwp, izp, rho)
            %SPEC_QUAD Construct an instance of this class
            %   Detailed explanation goes here
            if nargin > 3
                rho = 1;
            end            
            obj@spec_interface(iwp, izp, rho);
            obj.gain = GAIN;
        end
        
        function M = supply(obj)
            %SUPPLY quadratic performance specification
            gamma = obj.gain;
            Mu = -eye(obj.nwp)*gamma;
            My = eye(obj.nzp)/(gamma);

            M = blkdiag(My, Mu);
        end

       function [obj] = set_p(obj, p)
            %SET_P set a parameter when performing bisection
            %
            %
            %Example: Peak-to-Peak norm certifier
            %or l2 Gain bound

            obj.gain = p;
            
        end
    end
end

