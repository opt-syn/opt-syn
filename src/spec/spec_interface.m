classdef spec_interface
    % SPEC_INTERFACE Base class for a performance specification of an algorithm.
    %
    % A specification encodes a desired input-output property of the
    % algorithmic interconnection :math:`(F, G)` as a quadratic supply rate
    % on the performance channel :math:`(w_p, z_p)`. Concrete specifications
    % (stability, :math:`\ell_2`-gain, passivity, ...) subclass this
    % interface and override :mat:meth:`supply` / :mat:meth:`supply_quad`.
    %
    % The default (this base class) is a trivial stability specification with
    % an empty supply rate. The supported specialisations include the
    % energy-to-energy (:math:`\ell_2`) gain, the energy-to-peak
    % (generalized :math:`H_2`) gain, the peak-to-peak gain, and covariance
    % amplification (:math:`H_2`).
    %
    % .. note::
    %
    %    None of these specifications involve loop transformations, which
    %    keeps the assembly routines simple.

    properties
        iwp = [];       % Indices of the performance-input channel :math:`w_p`.
        izp = [];       % Indices of the performance-output channel :math:`z_p`.
        vars;           % Decision variables introduced by the specification.
        id = 0;         % Identifier / index of the specification.
        target = false; % Whether this specification is the optimization target
                        % within bisection iterations.
        rho = 1;        % Exponential discount factor :math:`\rho` of the specification.
        type = 'generic'; % Specification type tag.
    end

    methods
        function obj = spec_interface(iwp, izp)
            % SPEC_INTERFACE Construct a specification on a performance channel.
            %
            % :param iwp: Performance-input indices in the network.
            % :type iwp: double (vector)
            % :param izp: Performance-output indices in the network.
            % :type izp: double (vector)
            % :returns: A new specification object.
            % :rtype: spec_interface

            if nargin > 1
                obj.iwp = iwp;
                obj.izp = izp;
            end

        end


        function nzzp = nzp(obj)
            % NZP Number of performance outputs.
            %
            % :returns: Length of the :math:`z_p` index vector.
            % :rtype: int
            nzzp = length(obj.izp);
        end

        function nwwp = nwp(obj)
            % NWP Number of performance inputs.
            %
            % :returns: Length of the :math:`w_p` index vector.
            % :rtype: int
            nwwp = length(obj.iwp);
        end

        function [vars_spec, cons] = create_vars(obj, cons, name, config)
            % CREATE_VARS Create the decision variables for the specification.
            %
            % The base implementation introduces no variables. Subclasses
            % override this to add gain or multiplier variables and append
            % their defining constraints.
            %
            % :param cons: Accumulated LMI constraints.
            % :param name: Name suffix for the created variables.
            % :type name: char
            % :param config: Configuration options.
            % :type config: opt_config
            % :returns: ``[vars_spec, cons]`` — specification variables and the
            %    updated constraint set.
            % :rtype: struct, cell

            vars_spec = [];
        end

        function [M] = supply(obj, vars_spec)
            % SUPPLY Quadratic supply-rate matrix of the specification.
            %
            % Returns the matrix :math:`M` defining the running quadratic
            % cost on :math:`[z_p; w_p]`. The base implementation returns an
            % empty matrix (no constraint).
            %
            % :param vars_spec: Specification variables.
            % :type vars_spec: struct
            % :returns: Quadratic running-cost matrix :math:`M`.
            % :rtype: double
            M = [];
        end



        function [quad, objective] = supply_quad(obj, vars_spec)
            % SUPPLY_QUAD Decomposed quadratic supply rate.
            %
            % Splits the supply-rate matrix returned by :mat:meth:`supply`
            % into the output/input blocks used to assemble the LMI. If the
            % supply rate is empty, an empty :math:`\texttt{quad\_param}` is
            % returned.
            %
            % :param vars_spec: Specification variables.
            % :type vars_spec: struct
            % :returns: ``[quad, objective]`` — the decomposed quadratic
            %    (``quad_param``) and the objective contribution (``0`` for the
            %    base class).
            % :rtype: quad_param, double

            M = obj.supply(vars_spec);
            objective = 0;
            if isempty(M)
                quad = quad_param();                
            else
                nzp = length(obj.izp);
                nwp = length(obj.iwp);
                quad = quad_objective_decomp(M, 1:nzp, nzp + (1:nwp));
                
            end
            
        end
        

        function [obj] = set_p(obj, p)
            % SET_P Set the bisection parameter of the specification.
            %
            % Used by the bisection routine to update the swept quantity
            % (for example the convergence rate :math:`\rho` or an
            % :math:`\ell_2`-gain bound). The base implementation is a no-op.
            %
            % :param p: New value of the swept parameter.
            % :type p: double
            % :returns: The updated specification.
            % :rtype: spec_interface

        end

    end

    methods (Abstract)
        % create_vars(obj, cons)
        % supply(obj)
    end
        
        % function [vars, cons, iqc] = create_iqc(obj, cons)
        %     %CREATE_IQC Summary of this method goes here
        %     %   Detailed explanation goes here
        % 
        %     if nargin < 2;
        %         cons = [];
        %     end
        %     nwp = length(obj.iwp);
        %     nzp = length(obj.izp);
        %     switch obj.type
        %         case 'e2e'
        %             vars = []; iqc = iqc_e2e(nwp, nzp, obj.bound);
        %         case 'finite_l2'
        %             mu_l2 = lmim('mu_l2', 1, 1);
        % 
        %             cons = append_lmi(cons, mu_l2, obj.LMILAB);
        % 
        %             % cons = append_lmi(cons, mu_l2, obj.LMILAB);
        %             cons = append_lmi(cons, obj.finite_l2_bound - mu_l2, obj.LMILAB);
        % 
        %             vars = struct('mu_l2', mu_l2); 
        %             % iqc = iqc_finite_l2(nwp, nzp, mu_l2);
        %             % iqc = iqc_finite_l2(nwp, nzp, 2500);
        %             % iqc = iqc_e2e(nwp, nzp, 50);
        %             iqc = iqc_finite_l2(nwp, nzp, 50);
        %         otherwise
        %             vars = []; iqc = [];
        %     end
        % end
    
end
