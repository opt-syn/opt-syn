classdef spec_p2p < spec_interface
    %SPEC_P2P specification for a peak to peakgain
    %
    %sup_{T >= 0} norm(zp_k, 2)/norm(wp_k, 2) < gain
    %
    %
    % This is the peak-to-peak induced norm (*-norm of a linear system)



    properties
        type='p2p';
        gain = 1;
    end

    methods
        function obj = spec_p2p(GAIN, iwp, izp, rho)
            %SPEC_P2P Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 4
                rho = 1;
            end            
            obj@spec_interface(iwp, izp, rho);
            obj.gain = GAIN;
        end

        function M = supply(obj, vars_spec)
            %SUPPLY quadratic performance specification
            %
            %if the gain is fixed
            gamma = obj.gain;
            mu = vars_spec.mu_p2p;
            rrecip = obj.rho/(1-obj.rho);
            Mu = eye(obj.nwp) * rrecip * (gamma - mu);
            My = -eye(obj.nzp) * rrecip / gamma;

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

        function [vars, cons] = create_vars(obj, cons)
            %CREATE_VARS form the variables for the problem                        
            mu_p2p = lmim('mu_p2p', 1, 1);
            gam_p2p = lmim('gam_p2p', 1, 1);
            vars = struct('mu_p2p', mu_p2p, 'gam_p2p', gam_p2p);
            if ~obj.target
                cons = append_lmi(cons, obj.gain - gam_p2p, obj.LMILAB);
            end
            cons = append_lmi(cons, gam_p2p - mu_p2p, obj.LMILAB);
        end
    end
end

