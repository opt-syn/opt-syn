classdef spec_p2p < spec_interface
    % SPEC_P2P Peak-to-peak gain specification.
    %
    % Bounds the induced peak-to-peak norm (the :math:`*`-norm of a linear
    % system) from the performance input to the performance output: the peak
    % (over time) of the output is bounded by :math:`\gamma` times the peak
    % of the input,
    %
    % .. math::
    %
    %    \sup_{k \geq 0} \norm{z_{p,k}}_2 \; \leq \;
    %    \gamma \, \sup_{k \geq 0} \norm{w_{p,k}}_2 .

    properties        
        gain = 1;     % Peak-to-peak gain :math:`\gamma`.
        weight = 0.5; % Weighting between the peak and gain terms.
    end

    methods
        function obj = spec_p2p(GAIN, iwp, izp)
            % SPEC_P2P Construct a peak-to-peak gain specification.
            %
            % :param GAIN: Peak-to-peak gain bound :math:`\gamma`.
            % :type GAIN: double
            % :param iwp: Performance-input indices in the network.
            % :type iwp: double (vector)
            % :param izp: Performance-output indices in the network.
            % :type izp: double (vector)
            % :returns: A new peak-to-peak specification.
            % :rtype: spec_p2p
            % if nargin < 4
            %     rho = 1;
            % end            
            obj@spec_interface(iwp, izp);
            obj.gain = GAIN;
            obj.type = 'p2p';
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

        function quad = quad_terminal(obj, vars_spec)
            % QUAD_TERMINAL tbd
            %
            % :param vars_spec: Specification variables
            %    (expects ``mu_p2p``, ``gam_p2p``).
            % :type vars_spec: struct
            % :returns: Decomposed terminal quadratic.
            % :rtype: quad_param

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
            % SUPPLY_QUAD Decomposed quadratic supply rate.
            %
            % Builds the running quadratic on the performance input channel.
            % When this specification is the optimization target, the peak
            % gain ``gam_p2p`` is minimized.
            %
            % :param vars_spec: Specification variables
            %    (expects ``mu_p2p``, ``gam_p2p``).
            % :type vars_spec: struct
            % :returns: ``[quad, objective]`` — the decomposed quadratic and the
            %    objective contribution.
            % :rtype: quad_param, double

           
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
            % CREATE_VARS Create the peak-to-peak variables and constraints.
            %
            % Introduces ``mu_p2p`` and (when targeted) ``gam_p2p``, with the
            % ordering and non-negativity constraints of the peak-to-peak
            % certificate.
            %
            % :param cons: Accumulated LMI constraints.
            % :param name: Name suffix for the created variables.
            % :type name: char
            % :param config: Configuration options.
            % :type config: opt_config
            % :returns: ``[vars, cons]`` — struct with fields ``mu_p2p``,
            %    ``gam_p2p`` and the updated constraint set.
            % :rtype: struct, cell

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
