classdef spec_p2p < spec_interface
    %SPEC_P2P specification for a peak to peak gain
    %
    %In development

    %sup_{T >= 0} norm(zp_k, 2)/norm(wp_k, 2) < gain
    %
    %
    % This is the peak-to-peak induced norm (*-norm of a linear system)



    properties        
        gain = 1; %peak to peak gain
        weight = 0.5; %
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
            obj.type = 'p2p';
        end

        function M = supply(obj, vars_spec)
            %SUPPLY quadratic performance specification
            %
            %Args:
            %   vars_spec: problem variables in specification
            %
            %Returns:
            %   M: quadratic running cost matrix in the specification

            M0 = -obj.gain;

            M = kron(M0, eye(obj.nwp));
        end

        function quad = quad_terminal(obj, vars_spec)
            %QUAD_TERMINAL matrix for  the terminal p2p expression
            %
            %Args:
            %   vars_spec: problem variables in specification
            %
            %Returns:
            %   quad (quad_param): decomposed terminal quadratic specification

            %
            nwp = length(obj.iwp);
            nzp = length(obj.izp);

                
                Q0 =  -(vars_spec.mu_p2p - vars_spec.gam_p2p);
                U0 = vars_spec.gam_p2p;
                Q_term = drep(Q0, nwp);
                U_term = drep(U0, nzp);
 

            T_term = eye(nzp);
            S_term = zeros(nzp, nwp);

            
            quad = struct('Q', Q_term, 'U', U_term, 'T', T_term, 'S', S_term);
       end

       function [quad, objective] = supply_quad(obj, vars_spec)
            %SUPPLY_QUAD decomposed quadratic performance specification
            %
            %Args:
            %   vars_spec: problem variables in specification
            %
            %Returns:
            %   quad (quad_param): decomposed quadratic specification

           
               nwp = length(obj.iwp);
               
               Q_p2p = -drep(vars_spec.mu_p2p, nwp);
               
               S_p2p = zeros(nwp, 0);
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
            %Args:                        
            %   cons:  accumulated constraints            
            %   name: name of the specification
            %   config (opt_config): configuration options
            %
            %Returns:
            %   vars:  problem variables in specification
            %   cons:  accumulated constraints

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

