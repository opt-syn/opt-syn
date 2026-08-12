classdef spec_e2e < spec_interface
    % SPEC_E2E Energy-to-energy (:math:`\ell_2`) gain specification.
    %
    % Bounds the induced :math:`\ell_2`-gain from the performance input to
    % the performance output:
    %
    % .. math::
    %
    %    \limsup_{T \to \infty}
    %    \frac{\sum_{k=0}^{T} \norm{z_{p,k}}_2}{\sum_{k=0}^{T} \norm{w_{p,k}}_2}
    %    \; < \; \gamma .
    %
    % for any :math:`w_p` with finite energy. For a linear system this is 
    % the :math:`H_\infty` gain. It is the specification minimized when 
    % sweeping for the best achievable gain.

    properties        
        gain = 1; % :math:`\ell_2` gain :math:`\gamma` (hard constraint).
    end
    
    methods
        function obj = spec_e2e(GAIN, iwp, izp)
            % SPEC_E2E Construct an energy-to-energy gain specification.
            %
            % :param GAIN: Gain bound :math:`\gamma`.
            % :type GAIN: double
            % :param iwp: Performance-input indices in the network.
            % :type iwp: double (vector)
            % :param izp: Performance-output indices in the network.
            % :type izp: double (vector)
            % :returns: A new energy-to-energy gain specification.
            % :rtype: spec_e2e
   
            obj@spec_interface(iwp, izp);
            obj.gain = GAIN;
            obj.type='e2e';
        end
        
        function M = supply(obj, vars_spec)
            % SUPPLY Quadratic supply-rate matrix of the specification.
            %
            % Builds the block-diagonal supply rate
            % :math:`\operatorname{blkdiag}(\gamma^{-1} I,\, -\gamma I)` on
            % :math:`[z_p; w_p]`.
            %
            % :param vars_spec: Specification variables.
            % :type vars_spec: struct
            % :returns: Quadratic running-cost matrix :math:`M`.
            % :rtype: double

            gamma = obj.gain;

            Mzp = eye(obj.nwp)/gamma;            
            Mwp = -eye(obj.nzp)*(gamma);


            M = blkdiag(Mzp, Mwp);
        end

       function [quad, objective] = supply_quad(obj, vars_spec)
            % SUPPLY_QUAD Decomposed quadratic supply rate.
            %
            % When this specification is the optimization target, the gain
            % variable ``mu_l2`` is minimized; otherwise the base-class
            % decomposition is used.
            %
            % :param vars_spec: Specification variables (expects ``mu_l2``).
            % :type vars_spec: struct
            % :returns: ``[quad, objective]`` — the decomposed quadratic and the
            %    objective contribution.
            % :rtype: quad_param, double

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
            % SET_P Set the bisection parameter (the :math:`\ell_2` gain).
            %
            % :param p: New value of the gain bound.
            % :type p: double
            % :returns: The updated specification.
            % :rtype: spec_e2e

            obj.gain = p;
            
       end

       function [vars, cons] = create_vars(obj, cons, name, config)
            % CREATE_VARS Create the gain variable and its constraints.
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
            end
            % else
                % cons = append_lmi(cons, 1000 - mu_l2, config.LMILAB);
                cons = append_lmi(cons, mu_l2, config.LMILAB);
            
        end
    end
end
