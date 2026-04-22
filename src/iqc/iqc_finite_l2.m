classdef iqc_finite_l2 < iqc_loop_split
    %IQC_e2e Energy to Energy: Hinfinity gain for a linear system
    %
    % l2 -> l2 norm    
    %
    % norm(y)^2/norm(u)^2 <= gamma^2

    properties
        mu = 1;
    end
    
    methods
        function obj = iqc_finite_l2(nu, ny, mu)
            %IQC_HINF Construct an instance of this class
            %   Detailed explanation goes here

            if nargin < 3
                mu = 1;
            end

           
            %TODO: check signs
            Mu = -eye(nu)*mu;
            My = eye(ny);

            M = blkdiag(My, Mu);
            
            loop = [zeros(nu, ny), eye(ny);
                eye(nu), zeros(ny, nu)];

            obj@iqc_loop_split(eye(ny), M, loop, eye(nu), []);
            obj.mu = mu;
        end       
    end
end

