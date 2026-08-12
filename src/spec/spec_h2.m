classdef spec_h2 < spec_interface
    % SPEC_H2 Stochastic stability specification with input/output indices.
    %
    % The performance input $w_p$ is an i.d.d. zero-mean stochastic process
    % with bounded covariance $\mathbb{E}[||w_{p,k}||_2^2] \leq \sigma^2$
    % Imposes a stochastic stability constraint
    %
    % .. math::
    %
    %    \limsup_{T \rightarrow \infty} \frac{1}{T} \sum_{k=0}^{T}
    %    \mathbb{E}[|| z_{p,k}||^2_2] \leq \gamma^2 \sigma^2
    %
    % This is analogous to a primal H2 norm in the linear system setting.
    %
    %
    % :note: The h2 implementation requires that the stochastic input and
    % the oracle output are conditionally independent at each time k (zero
    % feedthrough in $D_{z w_p}$).
    %

    properties        
        gain = 10; % stochastic gain
        cov = 1; % covariance bound (scalar or matrix)
    end

    % type='passivity';
    % search_type = 'in';
    
    methods
        function obj = spec_h2(gain, cov, iwp, izp)
            % SPEC_PASSIVITY Construct a passivity specification.
            %
            % :param GAIN: Gain bound :math:`\gamma`.
            % :type GAIN: double
            % :param cov: Covariance of noise process
            % :type cov: double or symmetric matrix of double
            % :param iwp: Performance-input indices in the network.
            % :type iwp: double (vector)
            % :param izp: Performance-output indices in the network.
            % :type izp: double (vector)
            % :returns: A new passivity specification.
            % :rtype: spec_passivity
            % :raises: Error if ``iwp`` and ``izp`` have different lengths.
          
            obj@spec_interface(iwp, izp);
            obj.gain = gain;
            obj.cov= cov;
            obj.type = 'h2';
 
        end
        
        function M = supply(obj, vars_spec)
            % SUPPLY Quadratic supply-rate matrix for the passivity indices.
            %
            % :param vars_spec: Specification variables.
            % :type vars_spec: struct
            % :returns: Quadratic running-cost matrix :math:`M`.
            % :rtype: double
            
            M = -eye(obj.nwp);
        end

       function [quad, objective] = supply_quad(obj, vars_spec)
            % SUPPLY_QUAD Decomposed quadratic supply rate for h2 synthesis
            %
            %
            % :param vars_spec: Specification variables (expects ``ind_pass``).
            % :type vars_spec: struct
            % :returns: ``[quad, objective]`` — the decomposed quadratic and the
            %    objective contribution.
            % :rtype: quad_param, double


            % Q = -eye(obj.nwp);
            % S = zeros(0, obj.nwp);
            Q = [];
            S = zeros(obj.nzp,0);
            U = eye(obj.nzp);
            T = eye(obj.nzp);
            quad = quad_param(Q, S, U, T);

           if obj.target
  
              objective  = trace(vars_spec.Z);

           else
               objective = 0;

           end
       end

       function Omega = get_cov(obj)
           %GET_COV get the covariance matrix

           if isscalar(obj.cov)
               Omega = obj.cov * eye(obj.nzp);
           else
               Omega = obj.cov;
           end 
       end

       function objective = get_objective(obj, vars)
           %get the objective in the h2 problem
           if obj.target
               objective  = trace(vars.Z);

           else

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


            obj.gain = p;
            
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

            nzp = obj.nzp;
            if config.LMILAB                
                Z= lmim(['Z', name], nzp, nzp);            
            else
                Z= sdpvar(nzp, nzp);            
            end
            
            vars = struct('Z', Z);    

            cons = append_lmi(cons, Z, config.LMILAB);
            % if ~obj.target
                cons = append_lmi(cons, obj.gain^2 - trace(Z), config.LMILAB);
            % end
        end
    end
end
