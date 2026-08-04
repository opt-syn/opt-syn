classdef spec_quad < spec_interface
    % SPEC_QUAD General quadratic performance specification.
    %
    % Imposes a single fixed quadratic supply rate on the performance
    % channel :math:`(w_p, z_p)`:
    %
    % .. math::
    %
    %    \sum_{k=0}^{T} \begin{pmatrix} z_{p,k} \\ w_{p,k} \end{pmatrix}^\top
    %    \begin{pmatrix} Q & S \\ S^\top & R \end{pmatrix}
    %    \begin{pmatrix} z_{p,k} \\ w_{p,k} \end{pmatrix} \preceq 0 .
    %
    % The named specifications (:mat:class:`spec_l2`, :mat:class:`spec_e2e`,
    % :mat:class:`spec_passivity`, ...) are special cases obtained by fixing
    % the multiplier :math:`\left(\begin{smallmatrix} Q & S \\ S^\top & R \end{smallmatrix}\right)`.
    % Use this class directly when a desired supply rate is not covered by
    % one of the named specifications.

    properties
        M; % Quadratic performance matrix on :math:`[z_p; w_p]`.
    end
    
    methods
        function obj = spec_quad(M, iwp, izp)
            % SPEC_QUAD Construct a quadratic performance specification.
            %
            % :param M: Quadratic supply-rate matrix on :math:`[z_p; w_p]`.
            % :type M: double
            % :param iwp: Performance-input indices in the network.
            % :type iwp: double (vector)
            % :param izp: Performance-output indices in the network.
            % :type izp: double (vector)
            % :returns: A new quadratic specification.
            % :rtype: spec_quad
            obj@spec_interface(iwp, izp);
            obj.M = M;
            obj.type = 'quad';
        end
        
        function [M] = supply(obj)
            % SUPPLY Quadratic supply-rate matrix of the specification.
            %
            % Returns the stored matrix ``M`` unchanged.
            %
            % :returns: Quadratic running-cost matrix :math:`M`.
            % :rtype: double

            M = obj.M;
        end
    end
end
