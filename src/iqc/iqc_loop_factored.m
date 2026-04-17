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

            % psi_out = blkdiag(obj.Psi1, eye(nu2));

            Af = blkdiag(obj.Psi1.A, obj.Psi2.A);
            Bf = blkdiag(obj.Psi1.B, obj.Psi2.B);            

            Cf1 = [obj.Psi1.C, obj.C3];
            Df1 = [obj.Psi1.D, obj.D3];

            Cf2 = zeros(nu2, nx1 + nx2);
            Df2 = [zeros(nu2, nu1), eye(nu2)];

            Cf = [Cf1; Cf2];
            Df = [Df1; Df2];

            psi_out = ss(Af, Bf, Cf, Df, 1);
            
        end

        function psi_inv = psi2_inv(obj)
            %PSI2_INV inverse of the parameter psi2
            psi_inv = inv(obj.Psi2);
            %TODO: figure out clever conjugation
            %
            %how this should play with the zeros
        end


        function iqc_lift = lift(obj, c)
            %LIFT: kronecker by c (coordinates)


            Psi1_lift = ss_kron_eye(obj.Psi1, c);
            Psi2_lift = ss_kron_eye(obj.Psi2, c);
            C3_lift = kron_eye(obj.C3, c);
            D3_lift = kron_eye(obj.D3, c);
            loop_lift = kron(obj.loop, eye(c));
            M_lift = kron_eye(obj.M, c);
            X_lift = kron_eye(obj.X, c);
            
            iqc_lift = iqc_loop_split(Psi1_lift, M_lift, C3_lift, D3_lift, Psi2_lift, X_lift, loop_lift);
        end
    end
end

