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
            %NF number of states in the IQC filter
            nf_out = obj.nx1 + obj.nx2;            
        end

        function nx1_out = nx1(obj)
            %NX1 number of states in filter 1
            nx1_out = length(obj.Psi1.A);
        end

        function nx2_out = nx2(obj)
            %NX2 number of states in filter 2
            if ~isscalar(obj.Psi2)
                nx2_out = length(obj.Psi2.A);
            else
                nx2_out = 0;
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

        %% factorization routines
        function iqc_factored = factor(obj, perturb)
            %FACTOR spectral factorization of the IQC for synthesis
            %
            %

            if nargin < 2
                perturb = 1e-4;
            end
            
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

                iqc_factored = obj.perform_factorization(perturb );
                              
            else
                %the filter is already factored. fill in the information.
                C3 = zeros(obj.nx2, obj.nz);
                D3 = zeros(obj.nw, obj.nz);
                iqc_factored = iqc_loop_factored(Psi1, Psi2,...
                    C3, D3, obj.M, obj.X, obj.loop);
            end
        end

        function iqc_factored = perform_factorization(obj, perturb);
            %PERFORM_FACTORIZATION
            %the inner parts of the spectral factorization routines
            %
            % All factorization routines are based on the code of 
            % Lukas Schwenkel https://github.com/Schwenkel/mpc-iqc
            %
            % in the paper
            %
            %@article{Schwenkel2025,
            %   title={Output-feedback model predictive control under dynamic uncertainties using integral quadratic constraints},
            %   author={L. Schwenkel and J. K{\"o}hler and M. A. M{\"u}ller and F. Allg{\"o}wer}},
            %   year={2025},
            %   journal={arxiv:2504.00196},
            %   doi = {10.48550/arXiv.2504.00196},
            % }
            %
            % full credit to the authors

            %first determine the type of factorization
            
            if nargin < 2
                perturb = 1e-4;
            end                                
            np = obj.np;

            %top corner: select passive factorization
            M11 = obj.M(1:np, 1:np);
            

            if all(M11 == 0, "all")
                %passive factorization
                % iqc_factored = obj.passive_factorization();

                iqc_factored = obj.passive_factorize();
            else
                %PN-factorization
                %check the inertia
                % eM = eig(M);

                % epos= sum(eM>0);
                % eneg = sum(eM<0);
                % 
                % M_perturb = M;
                % if epos < nq
                %     M_perturb(1:nq, 1:nq) = M_perturb(1:nq, 1:nq) + eye(nq)*perturb;
                % end
                % 
                % obj.M = M_perturb 
                iqc_factored = obj.hinf_factorize();
            end

            %now do the factorization

            
        end

        function [iqc_factored, Vh, Z] = passive_factorize(obj)
            %PASSIVE_FACTORIZE perform an passive-infinity type factorization

            % 1. Constructing Psih_1
            % nq = obj.nq;
            % np = obj.np;

            %index into the relevant system            
            Psi1 = [obj.Psi1; zeros(obj.nq, obj.nz)];
            Psi2 = [ zeros(obj.np, obj.nw); obj.Psi2];

            nq = size(Psi1.B,2);
            np = size(Psi2.B,2);

            %form an equivalent supply rate
            % M11 = obj.M(1:np, 1:np);

            %form a generic supply rate for the Psi1 term
            M =  obj.M;
            MI = [eye(obj.np), zeros(obj.np, obj.nq); zeros(obj.nq, obj.np), zeros(obj.nq)];

            Q = Psi1.C'*MI*Psi1.C;
            R = Psi1.D'*MI*Psi1.D;
            S = Psi1.C'*MI*Psi1.D;

            %get a certificate
            [Zu, Ku, n0] = dare_anti(Psi1.A,Psi1.B,Q,R,S);


            D11hatu = chol(Psi1.B'*Zu*Psi1.B+R);
            C11hatu = D11hatu'\(Psi1.B'*Zu*Psi1.A+S');
            
            %add extra poles at zero to compensate for the lack of inverse
            %in discrete time
            Psi11u = ss(Psi1.A,Psi1.B,C11hatu,D11hatu,-1);                       
            
            Psi11 = tf('z')^(-n0)*Psi11u; 

            %the product (Psi11' Psi1) in state space
            Psi11Psi1 = minreal([Psi11; Psi1],[],false);
            Psi11Psi1 = balreal(Psi11Psi1);

            %matrices  for the product Psi1
            A1hat = Psi11Psi1.A;              B1hat = Psi11Psi1.B; 
            C11hat = Psi11Psi1.C(1:nq,:);     D11hat = Psi11Psi1.D(1:nq,:);   
            C1hat = Psi11Psi1.C(nq+1:end,:);

            %FLAG 1: the first transfer system
            Psi1h = ss(A1hat, B1hat, C11hat, D11hat, 1);



            % 2. Constructing Psih_12
            B1inv = B1hat/D11hat;
            D1inv = Psi1.D/D11hat;

            Psi1Psi11inv = ss(A1hat-B1inv*C11hat,B1inv,C1hat-D1inv*C11hat,D1inv,-1);
            [~, Psi11Psi11invMPsi2] = isproper(Psi1Psi11inv'*obj.M*Psi2);

            Psi11Psi11invMPsi2 = dss2ss(Psi11Psi11invMPsi2);
            
            Psi12Psi2 = minreal([Psi11Psi11invMPsi2; Psi2],[],false);
            Psi12Psi2 = balreal(Psi12Psi2);

            %data associated with Psih12
            A2hat = Psi12Psi2.A;              B2hat = Psi12Psi2.B; 
            C12hat = Psi12Psi2.C(1:nq,:);     D12hat = Psi12Psi2.D(1:nq,:);   
            C2hat = Psi12Psi2.C(nq+1:end,:);


            C3 = C12hat;
            D3 = D12hat;



            % 3. Constructing Psih_22
            %equivalent supply
            Q2 = Psi12Psi2.C'*blkdiag(eye(nq),-obj.M)*Psi12Psi2.C;
            R2 = Psi12Psi2.D'*blkdiag(eye(nq),-obj.M)*Psi12Psi2.D;
            S2 = Psi12Psi2.C'*blkdiag(eye(nq),-obj.M)*Psi12Psi2.D;
            Q2 = (Q2+Q2')/2;
            R2 = (R2+R2')/2;

            %stabilizing solution to riccati
            [Zs, Ks, Ls] = idare(A2hat,B2hat,Q2,R2,S2);
            D22hat = chol(B2hat'*Zs*B2hat+R2);
            C22hat = D22hat'\(B2hat'*Zs*A2hat+S2');

            %Flag Psi2
            Psi2h = ss(A2hat, B2hat, C22hat, D22hat, 1);
            
            % 4. package up the IQC
            Ahat = blkdiag(A1hat,A2hat);
            Bhat = blkdiag(B1hat,B2hat);
            Chat = [C11hat                   C12hat;
                    zeros(np,length(A1hat))  C22hat];

            Psih = struct;
            Psih.A = Ahat;
            Psih.B = Bhat;
            Psih.Chat = [C1hat C2hat];


    

            % 5. Computing Z, Xh, and Vh

            %get the state transformation/compression
            Vh = obj.compute_Vhat(Psih);
            Xh_V = Vh'*obj.X*Vh;

            Mhat = blkdiag(eye(nq),-eye(np));
            Q = [Chat; C1hat C2hat]'*blkdiag(Mhat,-M)*[Chat; C1hat C2hat];
            Q = (Q+Q')/2; % ensure symmetry
            Z = dlyap(Ahat',Q);
            


            Xhat = Xh_V+Z;
            Xhat = (Xhat+Xhat')/2;


            %package it all up
            iqc_factored = iqc_loop_factored(Psi1h, Psi2h,...
                    C3, D3, Mhat, Xhat, obj.loop);
        end
        

        function [iqc_factored, Vh, Z] = hinf_factorize(obj)
            %HINF_FACTORIZE perform an h-infinity type factorization

            % 1. Constructing Psih_1
            % nq = obj.nq;
            % np = obj.np;

            %index into the relevant system            
            Psi1 = [obj.Psi1; zeros(obj.nq, obj.nz)];
            Psi2 = [ zeros(obj.np, obj.nw); obj.Psi2];

            nq = size(Psi1.B,2);
            np = size(Psi2.B,2);

            %form an equivalent supply rate
            % M11 = obj.M(1:np, 1:np);

            %form a generic supply rate for the Psi1 term
            M =  obj.M;
            MI = [eye(np), zeros(np, nq); zeros(nq, np), zeros(nq)];

            Q = Psi1.C'*MI*Psi1.C;
            R = Psi1.D'*MI*Psi1.D;
            S = Psi1.C'*MI*Psi1.D;

            %get a certificate
            [Zu, Ku, n0] = dare_anti(Psi1.A,Psi1.B,Q,R,S);


            D11hatu = chol(Psi1.B'*Zu*Psi1.B+R);
            C11hatu = D11hatu'\(Psi1.B'*Zu*Psi1.A+S');
            
            %add extra poles at zero to compensate for the lack of inverse
            %in discrete time
            Psi11u = ss(Psi1.A,Psi1.B,C11hatu,D11hatu,-1);                       
            
            Psi11 = tf('z')^(-n0)*Psi11u; 

            %the product (Psi11' Psi1) in state space
            Psi11Psi1 = minreal([Psi11; Psi1],[],false);
            Psi11Psi1 = balreal(Psi11Psi1);

            %matrices  for the product Psi1
            A1hat = Psi11Psi1.A;              B1hat = Psi11Psi1.B; 
            C11hat = Psi11Psi1.C(1:nq,:);     D11hat = Psi11Psi1.D(1:nq,:);   
            C1hat = Psi11Psi1.C(nq+1:end,:);

            %FLAG 1: the first transfer system
            Psi1h = ss(A1hat, B1hat, C11hat, D11hat, 1);



            % 2. Constructing Psih_12
            B1inv = B1hat/D11hat;
            D1inv = Psi1.D/D11hat;

            Psi1Psi11inv = ss(A1hat-B1inv*C11hat,B1inv,C1hat-D1inv*C11hat,D1inv,-1);
            [~, Psi11Psi11invMPsi2] = isproper(Psi1Psi11inv'*obj.M*Psi2);

            Psi12Psi2 = minreal([Psi11Psi11invMPsi2; Psi2],[],false);
            Psi12Psi2 = balreal(Psi12Psi2);

            %data associated with Psih12
            A2hat = Psi12Psi2.A;              B2hat = Psi12Psi2.B; 
            C12hat = Psi12Psi2.C(1:nq,:);     D12hat = Psi12Psi2.D(1:nq,:);   
            C2hat = Psi12Psi2.C(nq+1:end,:);


            C3 = C12hat;
            D3 = D12hat;



            % 3. Constructing Psih_22
            %equivalent supply
            Q2 = Psi12Psi2.C'*blkdiag(eye(nq),-obj.M)*Psi12Psi2.C;
            R2 = Psi12Psi2.D'*blkdiag(eye(nq),-obj.M)*Psi12Psi2.D;
            S2 = Psi12Psi2.C'*blkdiag(eye(nq),-obj.M)*Psi12Psi2.D;
            Q2 = (Q2+Q2')/2;
            R2 = (R2+R2')/2;

            %stabilizing solution to riccati
            [Zs, Ks, Ls] = idare(A2hat,B2hat,Q2,R2,S2);
            D22hat = chol(B2hat'*Zs*B2hat+R2);
            C22hat = D22hat'\(B2hat'*Zs*A2hat+S2');

            %Flag Psi2
            Psi2h = ss(A2hat, B2hat, C22hat, D22hat, 1);
            
            % 4. package up the IQC
            Ahat = blkdiag(A1hat,A2hat);
            Bhat = blkdiag(B1hat,B2hat);
            Chat = [C11hat                   C12hat;
                    zeros(np,length(A1hat))  C22hat];

            Psih = struct;
            Psih.A = Ahat;
            Psih.B = Bhat;
            Psih.Chat = [C1hat C2hat];


    

            % 5. Computing Z, Xh, and Vh

            %get the state transformation/compression
            Vh = obj.compute_Vhat(Psih);
            Xh_V = Vh'*obj.X*Vh;

            Mhat = blkdiag(eye(nq),-eye(np));
            Q = [Chat; C1hat C2hat]'*blkdiag(Mhat,-M)*[Chat; C1hat C2hat];
            Q = (Q+Q')/2; % ensure symmetry
            Z = dlyap(Ahat',Q);
            


            Xhat = Xh_V+Z;
            Xhat = (Xhat+Xhat')/2;


            %package it all up
            iqc_factored = iqc_loop_factored(Psi1h, Psi2h,...
                    C3, D3, Mhat, Xhat, obj.loop);
        end
        
    
    
    
        function [ Vh] = compute_Vhat(obj, Psih)

                %
            % find Vhat such that Vhat*Ahat=Apsi*Vhat, Vhat*Bhat=Bpsi, and
            % Chat=Cpsi*Vhat by solving system of linear equations.
            %
            %Author: Lukas Schwenkel, 2025
        
        
            [Apsi, Bpsi, Cpsi, Dpsi] = ssdata(obj.get_psi);

            %work this through
            npsih = length(Psih.A);
            npsi = length(Apsi);

            %state transformation matrix
            Kvh = [ kron(eye(npsih),Apsi)-kron(Psih.A',eye(npsi));
                     kron(eye(npsih),Cpsi);
                     kron(Psih.B',eye(npsi))                        ];
            
            ansvh = [ zeros(npsi*npsih,1); Psih.Chat(:); Bpsi(:)       ];

            Vh = Kvh \ ansvh; 
                   
            Vh = reshape(Vh,[npsi, npsih]);
            
            % Alternative way to compute Vhat
            % [Ah,~,Ch,T1] = obsvf(Psih.A, Psih.B, Psih.Chat);
            % T2 = obsv(Ah(nW+1:end,nW+1:end), Ch(:,nW+1:end));
            % T3 = obsv(Psi.A, Psi.C);
            % Vhat = T3\T2*T1(nW+1:end,:);


        end
    end
end

