classdef iqc_loop_factored
    %IQC_LOOP_FACTORED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Psi1 = 1;   %primal filter (output of nonlinearity)
        Psi2 = 1;   %dual filter   (input of nonlinearity)        
        C3 = 0;     %crossover D matrix
        D3 = 0;     %crossover D matrix
        M = 0;      %running cost
        X = 0;      %terminal cost
        loop = 0;
    end
    
    methods
        function obj = iqc_loop_factored(Psi1, Psi2, C3, D3, M, X, loop)
            %IQC_LOOP_FACTORED Construct an instance of this class
            %   Detailed explanation goes here
            obj.Psi1 = Psi1;
            obj.Psi2 = Psi2;
            obj.C3 = C3;
            obj.D3 = D3;
            obj.M = M;
            obj.X = X;
            if nargin > 6
                obj.loop = loop;
            end

            if isscalar(obj.Psi1)
                obj.Psi1 = ss(obj.Psi1);
            end

            if isscalar(obj.Psi2)
                obj.Psi2 = ss(obj.Psi2);
            end
        end
        
        function hf = nf(obj)
            %NF number of states in the filter
            
            hf = length(obj.Psi1.A)+ length(obj.Psi2.A);
        end

        function  psi_out = get_psi(obj)
            %GET_PSI form the triangular multiplier

            psi_out = blkdiag(obj.Psi1, obj.Psi2);

            [ny1, nx1] = size(obj.Psi1.C);
            [~, nx2] = size(obj.Psi2.C);
            [~, nu1] = size(obj.Psi1.D);
            [~, nu2] = size(obj.Psi2.D);
            

            psi_out.C(1:ny1, nx1 + (1:nx2)) = obj.C3;
            psi_out.D(1:ny1, nu1 + (1:nu2)) = obj.D3;

            
        end

        function  psi_out = get_psi_13(obj)
            %GET_PSI form the triangular multiplier
            %only [psi1, psi3; 0, I]          

            [ny1, nx1] = size(obj.Psi1.C);   
            [~, nx2] = size(obj.Psi2.C);
            [~, nu1] = size(obj.Psi1.D);                               
            [~, nu2] = size(obj.Psi2.D);    

            psi_out = blkdiag(obj.Psi1, eye(nu2));

            

            psi_out.C(1:ny1, nx1 + (1:nx2)) = obj.C3;
            psi_out.D(1:ny1, nu1 + (1:nu2)) = obj.D3;

            % C13(1: ny1+1) ;

            %TODO: finish this
            psi_out = [];
            
        end

        function psi_inv = psi2_inv(obj)
            %PSI2_INV inverse of the parameter psi2
            psi_inv = inv(obj.Psi2);
        end
    end
end

