classdef spec_e2e_free < spec_e2e
    %SPEC_E2E_FREE specification for an energy to energy gain
    %
    %sum_k^T norm(zp_k)/norm(wp_k) < gain
    %
    %uses a Schur complement parameterization to allow for optimization
    
    properties
        type='e2e_schur';
        gain = 1;
    end
    
    methods
        function obj = spec_e2e_free(GAIN, iwp, izp, rho)
            %SPEC_QUAD Construct an instance of this class
            %   Detailed explanation goes here
            if nargin > 3
                rho = 1;
            end            
            obj@spec_e2e(GAIN, iwp, izp, rho);
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

