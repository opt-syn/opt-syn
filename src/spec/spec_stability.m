classdef spec_stability < spec_interface
    % SPEC_STABILITY Linear convergence / exponential stability specification.
    %
    % Enforces exponential (linear) convergence of the algorithmic
    % interconnection at rate :math:`\rho \in (0, 1)`: for every initial
    % condition there is a fixed point :math:`x^*(x_0)` with
    %
    % .. math::
    %
    %    \mav{c}{x_k - x^*(x_0) \\ w_k - w^*(x_0) \\ z_k - z^*(x_0)}_2
    %    \leq \gamma_0 \, \rho^{k} \norm{x_0 - x^*(x_0)}_2,
    %    \qquad \forall k \in \N .
    %
    % The rate is stored in the inherited ``rho`` property and enters the
    % LMI through the :math:`\rho`-transformation of the plant. This
    % specification carries no performance channel of its own (``iwp`` and
    % ``izp`` are empty); when :math:`\rho < 1` it :math:`\rho`-weights every
    % other specification on the same problem.

    methods
        function obj = spec_stability(rho)
            % SPEC_STABILITY Construct a linear-convergence specification.
            %
            % :param rho: Target convergence rate :math:`\rho \in (0, 1]`
            %    (default: ``1``, i.e. plain stability).
            % :type rho: double
            % :returns: A new stability specification.
            % :rtype: spec_stability
            if nargin == 0
                rho = 1;
            end
            obj@spec_interface([], []);
            obj.type = 'stability';     
            obj.rho = rho;
        end       

    end
end
