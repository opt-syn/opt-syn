classdef iqc_loop_split
    %IQC_LOOP_SPLIT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Psi1 = 1; %primal filter (output of nonlinearity)
        Psi2 = 1; %dual filter   (input of nonlinearity)
        loop = []; %loop transformation (in lft/star product form)
        M=0;    %running cost
        X=0;    %terminal cost
    end
    
    methods
        function obj = iqc_loop_split(Psi1, M, loop, Psi2, X)
            %IQC_LOOP_SPLIT Construct an instance of this class
            %   Detailed explanation goes here

            
            obj.Psi1 = Psi1;
            obj.M = M;

            if nargin > 2
                obj.loop = loop;
            else
                np = size(obj.Psi1.D, 1);
                nq = size(M, 1) - np;
                obj.loop = [zeros(np, np), eye(nq,np);
                        eye(np, nq), zeros(nq, nq)];
            end

            if nargin > 3
                obj.Psi2 = Psi2;
                obj.X = X;
            end
        end
        
        function nf_out = nf(obj)
            %number of states in the IQC filter
            nf_out = length(obj.Psi1.A);
            if ~isscalar(obj.Psi2)
                nf_out = nf_out + length(obj.Psi2.A);
            end
        end

        function sys = get_psi(obj)
            %get the multiplier by block diagonalization
            sys = blkdiag(obj.Psi1, obj.Psi2);
        end


        function iqc_lift = lift(obj, c)
            %LIFT: kronecker by c (coordinates)


            Psi1_lift = ss_kron_eye(obj.Psi1, c);
            Psi2_lift = ss_kron_eye(obj.Psi2, c);
            loop_lift = kron(obj.loop, eye(c));
            M_lift = kron_eye(obj.M, c);
            X_lift = kron_eye(obj.X, c);
            
            iqc_lift = iqc_loop_split(Psi1_lift, M_lift, loop_lift, Psi2_lift, X_lift);
        end
    end
end

