classdef spec_passivity < spec_interface
    % SPEC_PASSIVITY Passivity specification with input/output indices.
    %
    % Imposes a (possibly indexed) passivity supply rate on the performance
    % channel: for all horizons :math:`T` with :math:`x = 0`,
    %
    % .. math::
    %
    %    \sum_{k=0}^{T} z_{p,k}^\top w_{p,k}
    %    \; \geq \;
    %    \sum_{k=0}^{T} \big( \nu_w \norm{w_{p,k}}^2 + \nu_z \norm{z_{p,k}}^2 \big),
    %
    % where :math:`\nu_w = \texttt{ind_w}` is the input passivity index and
    % :math:`\nu_z = \texttt{ind_z}` is the output passivity index. The plain
    % passivity supply rate is recovered with both indices set to zero.

    properties        
        ind_w = 0; % Input passivity index :math:`\nu_w`.
        ind_z = 0; % Output passivity index :math:`\nu_z`.
    end

    % type='passivity';
    % search_type = 'in';
    
    methods
        function obj = spec_passivity(ind_w, ind_z, iwp, izp)
            % SPEC_PASSIVITY Construct a passivity specification.
            %
            % :param ind_w: Input passivity index :math:`\nu_w`.
            % :type ind_w: double
            % :param ind_z: Output passivity index :math:`\nu_z`.
            % :type ind_z: double
            % :param iwp: Performance-input indices in the network.
            % :type iwp: double (vector)
            % :param izp: Performance-output indices in the network.
            % :type izp: double (vector)
            % :returns: A new passivity specification.
            % :rtype: spec_passivity
            % :raises: Error if ``iwp`` and ``izp`` have different lengths.
          
            obj@spec_interface(iwp, izp);
            obj.ind_w = ind_w;
            obj.ind_z = ind_z;

            if length(iwp) ~= length(izp)
                error('Lengths of input and output vectors for passivity specification are not equal')
            end
        end
        
        function M = supply(obj, vars_spec)
            % SUPPLY Quadratic supply-rate matrix for the passivity indices.
            %
            % :param vars_spec: Specification variables.
            % :type vars_spec: struct
            % :returns: Quadratic running-cost matrix :math:`M`.
            % :rtype: double


            M0 = -[-obj.ind_z, 1; 1, -obj.ind_w];

            M = kron(M0, eye(obj.nwp));
        end

       function [quad, objective] = supply_quad(obj, vars_spec)
            % SUPPLY_QUAD Decomposed quadratic supply rate.
            %
            % When this specification is the optimization target, a passivity
            % margin ``ind_pass`` is introduced and maximized; otherwise the
            % base-class decomposition is used.
            %
            % :param vars_spec: Specification variables (expects ``ind_pass``).
            % :type vars_spec: struct
            % :returns: ``[quad, objective]`` — the decomposed quadratic and the
            %    objective contribution.
            % :rtype: quad_param, double

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
            % SET_P Set the bisection parameter (the input passivity index).
            %
            % :param p: New value of the input passivity index :math:`\nu_w`.
            % :type p: double
            % :returns: The updated specification.
            % :rtype: spec_passivity


            obj.ind_w = p;
            
       end

       function [vars, cons] = create_vars(obj, cons, name, config)
            % CREATE_VARS Create the passivity margin variable.
            %
            % :param cons: Accumulated LMI constraints.
            % :param name: Name suffix for the created variable.
            % :type name: char
            % :param config: Configuration options.
            % :type config: opt_config
            % :returns: ``[vars, cons]`` — struct with field ``ind_pass`` and the
            %    (unchanged) constraint set.
            % :rtype: struct, cell

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
