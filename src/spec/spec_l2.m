classdef spec_l2< spec_interface
    % SPEC_L2 :math:`\ell_2` stability specification.
    %
    % Certifies that the performance input has bounded influence on the
    % state, i.e. the algorithm is :math:`\ell_2`-stable. This is not an
    % :math:`\ell_2`-gain: there is no performance output (``izp = []``), only
    % a bounded quadratic penalty on the performance input :math:`w_p`.
    %
    % When the discount rate :math:`\rho \in (0, 1)`, this specification
    % implies Input-to-State Stability (ISS): there exist gains
    % :math:`\gamma_x, \gamma_w \geq 0` with
    %
    % .. math::
    %
    %    \norm{x_k - x^*(x_0)}_2^2 \leq
    %    \gamma_x \, \rho^k \norm{x_0 - x^*(x_0)}
    %    + \gamma_w \max_{t \in 0, \ldots, k} \norm{w_{p,t}}_2^2,
    %    \qquad \forall k \in \N .

    properties                        
        gain = 1e4; % :math:`\ell_2` bound (large default, effectively a slack).
    end
    % search_type = 'in';
    
    methods
        function obj = spec_l2(iwp, MU)
            % SPEC_L2 Construct an :math:`\ell_2` stability specification.
            %
            % :param iwp: Performance-input indices in the network.
            % :type iwp: double (vector)
            % :param MU: Optional :math:`\ell_2` bound (default: ``1e4``).
            % :type MU: double
            % :returns: A new :math:`\ell_2` stability specification.
            % :rtype: spec_l2
            izp = [];
            obj@spec_interface(iwp, izp);
            
            
            if nargin == 2
                obj.gain = MU;                
            end

            obj.type = 'l2_stability';        
        end
        
        function M = supply(obj, vars_spec)
            % SUPPLY Quadratic supply-rate matrix of the specification.
            %
            % Returns :math:`-\texttt{gain} \cdot I` on the performance input
            % channel.
            %
            % :param vars_spec: Specification variables.
            % :type vars_spec: struct
            % :returns: Quadratic running-cost matrix :math:`M`.
            % :rtype: double

            M0 = -obj.gain;

            M = kron(M0, eye(obj.nwp));
        end

       function [quad, objective] = supply_quad(obj, vars_spec)
            % SUPPLY_QUAD Decomposed quadratic supply rate.
            %
            % When this specification is the optimization target, the bound
            % ``mu_l2`` is minimized; otherwise the base-class decomposition
            % is used.
            %
            % :param vars_spec: Specification variables (expects ``mu_l2``).
            % :type vars_spec: struct
            % :returns: ``[quad, objective]`` — the decomposed quadratic and the
            %    objective contribution.
            % :rtype: quad_param, double

           if obj.target
               nwp = length(obj.iwp);

               Q_l2 = -vars_spec.mu_l2;
               
               S_l2 = -zeros(nwp, 0)';
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
            % SET_P Set the bisection parameter (the :math:`\ell_2` bound).
            %
            % :param p: New value of the swept parameter, e.g. a convergence
            %    rate :math:`\rho` or an :math:`\ell_2`-gain bound.
            % :type p: double
            % :returns: The updated specification.
            % :rtype: spec_l2
            
            obj.gain = p;
           
       end

       function [vars, cons] = create_vars(obj, cons, name, config)
            % CREATE_VARS Create the :math:`\ell_2` bound variable and constraints.
            %
            % :param cons: Accumulated LMI constraints.
            % :param name: Name suffix for the created variable.
            % :type name: char
            % :param config: Configuration options.
            % :type config: opt_config
            % :returns: ``[vars, cons]`` — struct with field ``mu_l2`` and the
            %    updated constraint set.
            % :rtype: struct, cell

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
