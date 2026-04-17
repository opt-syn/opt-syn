function [S_lift] = ss_kron_eye(S, c)
%SS_KRON kronecker with the identity matrix
%is a lift of the system
%with the right permutation structure

A_lift = kron_eye(S.A, c);
B_lift = kron_eye(S.B, c);
C_lift = kron_eye(S.C, c);
D_lift = kron_eye(S.D, c);

Ts = S.Ts;

if isa(S, 'ss')
    S_lift = ss(A_lift, B_lift, C_lift, D_lift, Ts);
else
    S_lift = sdpss(A_lift, B_lift, C_lift, D_lift, Ts);
end

end

