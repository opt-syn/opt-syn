classdef spec_e2e < spec_interface
    %SPEC_E2E specification for an energy to energy gain
    
    %limsup_{T -> inf} sum_{k=0}^T norm(zp_k, 2)/norm(wp_k, 2) < gain
    %
    %
    % For a linear system, this is the H infinity gain 

    
    
    properties        
        gain = 1; %l2 gain (hard constraint)
    end
    
    methods
        function obj = spec_e2e(GAIN, iwp, izp)
            %SPEC_E2E Constructor
            %
            %Args:            
            %   iwp:    performance inputs in the network    
            %   izp:    performance outputs in the network            
   
            obj@spec_interface(iwp, izp);
            obj.gain = GAIN;
            obj.type='e2e';
        end
        
        function M = supply(obj, vars_spec)
            %SUPPLY quadratic performance specification
            %
            %Args:
            %   vars_spec: problem variables in specification
            %
            %Returns:
            %   M: quadratic running cost matrix in the specification

            gamma = obj.gain;

            Mzp = eye(obj.nwp)/gamma;            
            Mwp = -eye(obj.nzp)*(gamma);


            M = blkdiag(Mzp, Mwp);
        end

       function [quad, objective] = supply_quad(obj, vars_spec)
            %SUPPLY_QUAD decomposed quadratic performance specification            %
            %Args:
            %   vars_spec: problem variables in specification
            %
            %Returns:
            %   quad (quad_param): decomposed quadratic specification

           if obj.target
               nwp = length(obj.iwp);
               nzp = length(obj.izp);
               objective  = vars_spec.mu_l2;

               Q_e2e = -drep(vars_spec.mu_l2, nwp);
               T_e2e = eye(nzp);
               S_e2e = zeros(nzp, nwp);
               U_e2e = drep(vars_spec.mu_l2, nzp);

               quad = struct('Q', Q_e2e, 'T', T_e2e, 'S', S_e2e, 'U', U_e2e);

               

           else
               quad = supply_quad@spec_interface(obj, vars_spec);
               objective = 0;
           end
       end

       function [obj] = set_p(obj, p)
            %SET_P set a parameter when performing bisection
            %
            %Args:
            %   p: new value of the parameter (l2 gain)

            obj.gain = p;
            
       end

       function [vars, cons] = create_vars(obj, cons, name, config)
            %CREATE_VARS form the variables for the problem                        
            %
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
                mu_l2 = lmim(['mu_l2', name], 1, 1);            
            else
                mu_l2 = sdpvar(1, 1);            
            end
            
            vars = struct('mu_l2', mu_l2);
            if ~obj.target
                cons = append_lmi(cons, obj.gain - mu_l2, config.LMILAB);
            else
                % cons = append_lmi(cons, 1000 - mu_l2, config.LMILAB);
                cons = append_lmi(cons, mu_l2, config.LMILAB);
            end
        end
    end
end

