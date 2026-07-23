classdef spec_quad < spec_interface
    %SPEC_QUAD quadratic performance
    % :math:`\sum_{i=0}^T \begin{pmatrix} w_k \\ z_k \end{pmatrix}^\top
    % \begin{pmatrix} Q & S \\ S^\top & T^\top U^{-1} T \end{pmatrix}
    % \begin{pmatrix} w_k \\ z_k \end{pmatrix} \leq -\sum_{i=0}^T \epsilon
    % ||w_k||_2^2` 

    
    %
    %[zp]' [Q  S]  [zp] << 0
    %[wp]  [S' R]  [wp] 

    
    properties
        M; %quadratic performance matrix        
    end
    
    methods
        function obj = spec_quad(M, iwp, izp)
            %SPEC_QUAD Constructor
            obj@spec_interface(iwp, izp);
            obj.M = M;
            obj.type = 'quad';
        end
        
        function [M] = supply(obj)
            %SUPPLY quadratic performance specification
            %
            %Args:
            %   vars_spec: problem variables in specification
            %
            %Returns:
            %   M: quadratic running cost matrix in the specification

            M = obj.M;
        end
    end
end

