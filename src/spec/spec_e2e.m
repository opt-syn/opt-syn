classdef spec_e2e < spec_interface
    %SPEC_E2E specification for an energy to energy gain
    %
    %limsup_{T -> inf} sum_{k=0}^T norm(zp_k, 2)/norm(wp_k, 2) < gain
    %
    %
    % For a linear system, this is the H infinity gain 

    
    
    properties
        type='e2e';
        gain = 1;
    end
    
    methods
        function obj = spec_e2e(GAIN, iwp, izp)
            %SPEC_E2E Construct an instance of this class
            %   Detailed explanation goes here
            % if nargin < 4
            %     rho = 1;
            % end            
            obj@spec_interface(iwp, izp);
            obj.gain = GAIN;
        end
        
        function M = supply(obj, vars_spec)
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

       function [vars, cons] = create_vars(obj, cons, name, config)
            %CREATE_VARS form the variables for the problem                        
            if nargin < 3
                name = [];
            end

            if config.LMILAB
                mu_l2 = lmim(['mu_l2', name], 1, 1);            
            else
                mu_l2 = sdpvar(1, 1);            
            end
            
            vars = struct('mu_l2', mu_l2);
            if ~obj.target
                cons = append_lmi(cons, obj.gain - mu_l2, config.LMILAB);
            end
        end
    end
end

