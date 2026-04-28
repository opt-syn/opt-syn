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
                obj.Psi2 = sdpss(Psi2);
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
            if isnumeric(obj.Psi2.D)                
                nw_out = size(obj.Psi2.D, 2);
            else
                nw_out = dim(obj.Psi2.D, 2);
            end        
        end

        function nz_out = nz(obj)
            %number of inputs (channel one)
            if isnumeric(obj.Psi1.D)                
                nz_out = size(obj.Psi1.D, 2);
            else
                nz_out = dim(obj.Psi1.D, 2);
            end
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
            %number of outputs (channel two)
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


        function iqc_rec = recover(obj, lmi_out)
            %RECOVER: recover the numerical values of an IQC from an LMI
            %solution

            M_rec = double(double(obj.M, lmi_out));
            X_rec = double(double(obj.X, lmi_out));
            loop_rec = obj.loop;

            Cf1 = double(double(obj.Psi1.C, lmi_out));
            Cf2 = double(double(obj.Psi2.C, lmi_out));
            Df1 = double(double(obj.Psi1.D, lmi_out));
            Df2 = double(double(obj.Psi2.D, lmi_out));

            Psi1_rec = ss(obj.Psi1.A, obj.Psi1.B, Cf1, Df1, 1);
            Psi2_rec = ss(obj.Psi2.A, obj.Psi2.B, Cf2, Df2, 1);

            iqc_rec = iqc_loop_split(Psi1_rec, M_rec, loop_rec, Psi2_rec, X_rec);

            iqc_rec.Psi1 = Psi1_rec;
            iqc_rec.Psi2 = Psi2_rec;
        end        

        %% factorization routine
        function iqc_factored = factor(obj)
            %FACTOR spectral factorization of the IQC for synthesis
            %
            %
            
            Psi1 = ss(obj.Psi1.A, obj.Psi1.B, obj.Psi1.C, obj.Psi1.D, 1);
            Psi2 = ss(obj.Psi2.A, obj.Psi2.B, obj.Psi2.C, obj.Psi2.D, 1);

            %Do we need to factor at all?
            need_to_factor = true;
            dim_check = (obj.nw == obj.nq) && (obj.nz == obj.np) ;
            
            
            if dim_check
                %both filters should be stable, and Psi2 should have a
                %stable inverse
                if ssize(Psi1.A, 1) == 0
                    e1 = [];
                else
                    e1 = eig(obj.Psi1);
                end
                if ssize(Psi2.A, 1) == 0
                    e2 = [];
                    e2inv = [];

                    invscal = (min(svd(Psi2.D)) < 1e-8);
                else
                    e2 = eig(Psi2);
                    P2inv = inv(Psi2);
                    e2inv = eig(P2inv);

                    invscal = false;
                end

                eall = [e1; e2; e2inv];
                %TODO: check the comparator: >= or >?
                stab_check = any(abs(eall) > 1) || invscal;

                need_to_factor = stab_check;
            end

            if need_to_factor
                %TODO: implement the spectral factorization
                error('TODO: spectral factorization not yet implemented')
                iqc_factored = [];
            else
                %the filter is already factored. fill in the information.
                C3 = zeros(obj.nf, obj.nz);
                D3 = zeros(obj.nw, obj.nz);
                iqc_factored = iqc_loop_factored(Psi1, Psi2,...
                    C3, D3, obj.M, obj.X, obj.loop);
            end
        end
    end
end

