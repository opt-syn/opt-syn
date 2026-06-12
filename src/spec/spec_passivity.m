classdef spec_passivity < spec_interface
    %SPEC_PASSIVITY specification for an passivity specification
    %
    %limsup_{T -> inf} sum_{k=0}^T (zp_k)' (wp_k) > 
    %                               sum_{k=0}^T ind_w |wp_k|^2 + ind_z |zp_k|^2
    %


    
    
    properties
        type='passivity';
%         gain = 1;
        ind_w = 0;
        ind_z = 0;
        search_type = 'in';
    end
    
    methods
        function obj = spec_passivity(ind_w, ind_z, iwp, izp)
            %SPEC_pass Construct an instance of this class
            %   Detailed explanation goes here
            % if nargin < 4
            %     rho = 1;
            % end            
            obj@spec_interface(iwp, izp);
            obj.ind_w = ind_w;
            obj.ind_z = ind_z;

            if length(iwp) ~= length(izp)
                error('Lengths of input and output vectors for passivity specification are not equal')
            end
        end
        
        function M = supply(obj, vars_spec)
            %SUPPLY quadratic performance specification
            %for passivity indices

            M0 = [-obj.ind_z, 1; 1, -obj.ind_w];

            M = kron(M0, eye(obj.nwp));
        end

       function [quad, objective] = supply_quad(obj, vars_spec)
           %SUPPLY_QUAD decomposed quadratic performance specification

           if obj.target
               nwp = length(obj.iwp);
               nzp = length(obj.izp);
               

               Q_pass = drep(vars_spec.ind_pass, nwp);
               
               S_pass = -eye(nzp, nwp);

               if obj.ind_z ~= 0
                   T_pass = eye(nzp);
                   U_pass = 1/obj.ind_z;
               else
                   T_pass = zeros(0, nzp);
                   U_pass = [];
               end

               quad = struct('Q', Q_pass, 'T', T_pass, 'S', S_pass, 'U', U_pass);

               objective  = -vars_spec.ind_pass;

           else
               quad = supply_quad@spec_interface(obj, vars_spec);
               objective = 0;

           end
       end

       function [obj] = set_p(obj, p)
            %SET_P set a parameter when performing bisection


            obj.ind_w = p;
            
       end

       function [vars, cons] = create_vars(obj, cons, name, config)
            %CREATE_VARS form the variables for the problem                        
            if nargin < 3
                name = [];
            end

            if config.LMILAB
                ind_pass = lmim(['ind_pass', name], 1, 1);            
            else
                ind_pass = sdpvar(1, 1);            
            end
            
            vars = struct('ind_pass', ind_pass);            
        end
    end
end

