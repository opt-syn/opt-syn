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
        
        function obj = factor(obj)
            %the IQC is already factored, do nothing
        end

        %% dimension counters

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

        function hf = nf(obj)
            %NF number of states in the filter            
            hf = length(obj.Psi1.A)+ length(obj.Psi2.A);
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
        
                    %care must be taken with empties
                    if isempty(b.C3)
                        if isempty(obj.C3)
                            C3 = zeros(n1+n2, 0);
                        else
                            C3 = [obj.C3, zeros(n2, size(obj.C3, 1))];
                        end
                    else
                        if isempty(obj.C3)
                            C3 = [zeros(n1, size(b.C3, 1)); b.C3, ];
                        else
                            C3 =  blkdiag(obj.C3, b.C3);
                        end
                    end          



                    D3 = blkdiag(obj.D3, b.D3);
    
                    iqc = iqc_loop_factored(Psi1, Psi2, C3, D3, M, X, loop);

                end
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

        function G_wrap = wrap_synth(obj, G, n)
            %wrap a genearlized plant for controller synthesis

            %
            % [z] = G [w]
            % [y]     [u]
            %to [Psi1, Psi3 Psi2^-1]  [G Psi2^-1]
            %   [0,               I]  [I]
            %   [
            %with G cl



            %routine copied and modified from 
            % https://github.com/Schwenkel/mpc-iqc/blob/main/multi-objective-iqc-synthesis/syn_step.m
            %by Lukas Schwenkel

            %dimension counts

            n1 = length(obj.Psi1.A);
            n2 = length(obj.Psi2.A);
            nx = length(G.A);    
            
            
            nw = n.nw;            
            nz = n.nz;
            
            nwp = n.nwp;
            nzp = n.nzp;
            
            nu = n.nu;            
            ny = n.ny;

            %index groups (unfortunately, the nature of the genplant is
            %lost. fix with overrides?)
            p = 1:nw;
            q = 1:nz;

            w = nw + (1:nwp);
            z = nz + (1:nzp);

            u = nw + nwp + (1:nu);
            y = nz + nzp + (1:ny);

            C2i = -obj.Psi2.D\obj.Psi2.C;
            A2i = obj.Psi2.A+obj.Psi2.B*C2i;
            B2i = obj.Psi2.B/obj.Psi2.D;
            
            
            %TODO: check the number of outputs here.
            A = [ obj.Psi1.A        obj.Psi1.B*G.D(q,p)*C2i   obj.Psi1.B*G.C(q,:)
                  zeros(n2,n1)   A2i                    zeros(n2,nx)
                  zeros(nx,n1)   G.B(:,p)*C2i           G.A          ];
            B = [ obj.Psi1.B*(G.D(q,p)/obj.Psi2.D)  obj.Psi1.B*G.D(q,w)  obj.Psi1.B*G.D(q,u)
                  B2i                         zeros(n2,nwp)      zeros(n2,nu)
                  G.B(:,p)/obj.Psi2.D            G.B(:,w)          G.B(:,u)   ];

            Dsp = (obj.Psi1.D*G.D(q,p)+obj.D3);

            C = [ obj.Psi1.C       obj.C3+Dsp*C2i  obj.Psi1.D*G.C(q,:)
                  zeros(nzp,n1)  G.D(z,p)*C2i     G.C(z,:);
                  zeros(ny,n1)  G.D(y,p)*C2i     G.C(y,:) ];
            D = [ Dsp/obj.Psi2.D       obj.Psi1.D*G.D(q,w)  obj.Psi1.D*G.D(q,u)
                  G.D(z,p)/obj.Psi2.D  G.D(z,w)          G.D(z,u)  
                  G.D(y,p)/obj.Psi2.D  G.D(y,w)          G.D(y,u) ];

            P = ss(A, B, C, D, 1);


            %dimensions don't change: square multiplier
            G_wrap = genplant(P, n);

        end
    end
end

