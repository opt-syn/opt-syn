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
        weight = 0.5; %used for the terminal cost (theorem 3 of Schwenkel 2026)
    end

    methods
        function obj = spec_p2p(GAIN, iwp, izp)
            %SPEC_P2P Construct an instance of this class
            %   Detailed explanation goes here
            % if nargin < 4
            %     rho = 1;
            % end            
            obj@spec_interface(iwp, izp);
            obj.gain = GAIN;
        end

        function M = supply(obj, vars_spec)
            %SUPPLY quadratic performance specification
            %for passivity indices

            M0 = -obj.gain;

            M = kron(M0, eye(obj.nwp));
        end

        function quad = quad_terminal(obj, vars_spec)
            %matrix for  the terminal p2p expression
            %
            alpha = obj.rho^2 / (1-obj.rho^2);
            nwp = length(obj.iwp);
            nzp = length(obj.izp);

            if obj.target
                Q0 = -obj.weight*alpha * (vars_spec.gam_p2p - vars_spec.mu_p2p);
                U0 = vars_spec.gam_p2p *  (1/(obj.weight*alpha));
                Q_term = drep(Q0, nwp);
                U_term = drep(U0, nzp);
            else
                
                Q0 = -obj.weight * alpha * (obj.gain - vars_spec.mu_p2p);
                U0 = obj.gain / (obj.weight*alpha);
                Q_term = drep(Q0, nwp);
                U_term = kron(U0, eye(nzp));
                
            end

            T_term = eye(nzp);
            S_term = zeros(nzp, nwp);

            
            quad = struct('Q', Q_term, 'U', U_term, 'T', T_term, 'S', S_term);
       end

       function [quad, objective] = supply_quad(obj, vars_spec)
           %SUPPLY_QUAD decomposed quadratic performance specification

           
               nwp = length(obj.iwp);
               
               Q_p2p = -vars_spec.mu_p2p;
               
               S_p2p = -zeros(nwp, 0);
               T_p2p = [];
               U_p2p = [];
  

               quad = struct('Q', Q_p2p, 'T', T_p2p, 'S', S_p2p, 'U', U_p2p);
        
           
           if obj.target
               objective  = vars_spec.gam_p2p;

           else
               
               objective = 0;

           end
       end

        function [vars, cons] = create_vars(obj, cons, name, config)
            %CREATE_VARS form the variables for the problem                        
            
            if nargin < 3
                name = [];
            end

            if config.LMILAB
                mu_p2p = lmim(['mu_p2p', name], 1, 1);
            else
                mu_p2p = sdpvar(1, 1);
            end
            
            
            if obj.target
                if config.LMILAB
                    gam_p2p = lmim(['gam_p2p', name], 1, 1);
                else
                    mu_p2p = sdpvar(1, 1);
                end
                cons = append_lmi(cons, obj.gain - gam_p2p, config.LMILAB);
            else
                gam_p2p = obj.gain;
            end

            vars = struct('mu_p2p', mu_p2p, 'gam_p2p', gam_p2p);
            cons = append_lmi(cons, gam_p2p - mu_p2p, config.LMILAB);
            cons = append_lmi(cons,  mu_p2p, config.LMILAB);
        end
    end
end

