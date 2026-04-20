classdef iqc_e2e < iqc_loop_split
    %IQC_e2e Energy to Energy: Hinfinity gain for a linear system
    %
    % l2 -> l2 norm    
    %
    % norm(y)^2/norm(u)^2 <= gamma^2

    properties
        gamma2
    end
    
    methods
        function obj = iqc_e2e(nu, ny, gamma2)
            %IQC_HINF Construct an instance of this class
            %   Detailed explanation goes here

            if nargin < 3
                gamma2 = 1;
            end

           
            %TODO: check signs
            Mu = -eye(nu)*gamma2;
            My = eye(ny);

            M = blkdiag(My, Mu);
            
            obj@iqc_loop_split(1, M, 1, 1, []);
            obj.gamma = gamma;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

