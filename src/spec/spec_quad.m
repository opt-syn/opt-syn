classdef spec_quad < spec_interface
    % SPEC_QUAD General quadratic performance specification.
    %
    %
    %
    %[zp]' [Q  S]  [zp] < 0
    %[wp]  [S' R]  [wp] 

    
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
