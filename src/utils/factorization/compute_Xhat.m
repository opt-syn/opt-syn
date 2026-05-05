function [Xh, Vh] = compute_Xhat(Psi1, Psi2, Psih, Z, X)
    %
    % find Vhat such that Vhat*Ahat=Apsi*Vhat, Vhat*Bhat=Bpsi, and
    % Chat=Cpsi*Vhat by solving system of linear equations.
    %
    %Author: Lukas Schwenkel, 2025

    
    npsih = length(Psih.A1)+length(Psih.A2);
    Apsi = blkdiag(Psi1.A, Psi2.A);
    Bpsi = blkdiag(Psi1.B, Psi2.B);
    Cpsi = [Psi1.C Psi2.C];
    npsi = length(Apsi);
    Vh = [ kron(eye(npsih),Apsi)-kron(Psih.A',eye(npsi));
             kron(eye(npsih),Cpsi);
             kron(Psih.B',eye(npsi))                        ] \ ...
           [ zeros(npsi*npsih,1); Psih.Chat(:); Bpsi(:)       ];
    Vh = reshape(Vh,[npsi, npsih]);
    
    % Alternative way to compute Vhat
    % [Ah,~,Ch,T1] = obsvf(Psih.A, Psih.B, Psih.Chat);
    % T2 = obsv(Ah(nW+1:end,nW+1:end), Ch(:,nW+1:end));
    % T3 = obsv(Psi.A, Psi.C);
    % Vhat = T3\T2*T1(nW+1:end,:);

    Xh = Vh'*X*Vh+Z;
    Xh = (Xh+Xh')/2;
end