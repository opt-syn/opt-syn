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

            
            obj.Psi1 = sdpss(Psi1);
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

            obj.Psi2 = sdpss(Psi2);
        end
        
        function nf_out = nf(obj)
            %number of states in the IQC filter
            nf_out = length(obj.Psi1.A);
            if ~isscalar(obj.Psi2)
                nf_out = nf_out + length(obj.Psi2.A);
            end
        end

        function nw_out = nw(obj)
            %number of inputs (channel two)
            nw_out = size(obj.Psi1.B, 2);            
        end

        function nz_out = nz(obj)
            %number of inputs (channel two)
            nz_out = size(obj.Psi2.B, 2);            
        end

        function np_out = np(obj)
            %number of outputs (channel one)
            if isnumeric(obj.Psi1.D)                
                np_out = size(obj.Psi1.D, 1);
            else
                np_out = dim(obj.Psi1.D, 1);
            end
        end

        function nq_out = nq(obj)
            %number of outputs (channel one)
            if isnumeric(obj.Psi2.D)
                nq_out = size(obj.Psi2.D, 1);
            else                
                nq_out = dim(obj.Psi2.D, 1);
            end
        end

        function sys = get_psi(obj)
            %get the multiplier by block diagonalization
            sys = blkdiag(obj.Psi1, obj.Psi2);
        end

        function iqc = blkdiag(obj, b)
            %BLKDIAG: block diagonal of the multipliers

            % if length(varargin) == 1
                % b = varargin{1};
                if isempty(b)
                    iqc = obj;
                else                   
                    Psi1 = blkdiag(obj.Psi1, b.Psi1);
                    Psi2 = blkdiag(obj.Psi2, b.Psi2);
                    X = blkdiag(obj.X, b.X);
        
        
                    n1 = obj.nw;
                    m1 = obj.nz;
        
                    n2 = b.nw;
                    m2 = b.nz;   
                    
    
                    [loop] = outer_blkdiag(obj.loop, b.loop, n1, m1, n2, m2);
    
                    % loop = E1*obj.loop*E1' + E2*b.loop*E2';
        
                    n1 = obj.np;
                    m1 = obj.nq;
                    n2 = b.np;
                    m2 = b.nq;
    
                    [M] = outer_blkdiag(obj.M, b.M, n1, m1, n2, m2);
        
    
                    iqc = iqc_loop_split(Psi1, M, loop, Psi2, X);

                end
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

