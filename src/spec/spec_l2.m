classdef spec_l2< spec_interface
    %SPEC_L2 specification for l2 stability
    %
    % Not L2 gain (no output zp)
    %
    % When rho in (0, 1), this implies Input to State Stability
    %


    
    
    properties
        type='passivity';        
        ind_w = 0;        
        search_type = 'in';
        gain = 1e4;
    end
    
    methods
        function obj = spec_l2(iwp, MU)
            %SPEC_l2 Construct an instance of this class
            %   Detailed explanation goes here
            % if nargin < 4
            %     rho = 1;
            % end            
            izp = [];
            obj@spec_interface(iwp, izp);
            
            
            if nargin == 2
                obj.gain = MU;                
            end
        end
        
        function M = supply(obj, vars_spec)
            %SUPPLY quadratic performance specification
            %for passivity indices

            M0 = -obj.gain;

            M = kron(M0, eye(obj.nwp));
        end

       function [quad, objective] = supply_quad(obj, vars_spec)
           %SUPPLY_QUAD decomposed quadratic performance specification

           if obj.target
               nwp = length(obj.iwp);
               
               

               Q_l2 = -vars_spec.mu_l2;
               
               S_l2 = -zeros(nwp, 0);
               % S_l2 = [];

               T_l2 = [];
               U_l2 = [];
  

               quad = struct('Q', Q_l2, 'T', T_l2, 'S', S_l2, 'U', U_l2);

               objective  = vars_spec.mu_l2;

           else
               quad = supply_quad@spec_interface(obj, vars_spec);
               objective = 0;

           end
       end

       function [obj] = set_p(obj, p)
            %SET_P set a parameter when performing bisection

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
            else
                % cons = append_lmi(cons, 1000 - mu_l2, config.LMILAB);
                cons = append_lmi(cons, mu_l2, config.LMILAB);
            end
        end
    end
end

