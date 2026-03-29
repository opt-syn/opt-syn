function [Pvar, psd_block] = positive_real_lmi(A2b, B2b, C2b, D2b, pass_tol, LMILAB)
%POSITIVE_REAL_LMI  enforce that a discrete-time linear system is
%positive-real

[n, m] = size(B2b);

if n == 0
    Pvar = [];
    if LMILAB
        p = dim(D2b, 1);
    else
        p = size(D2b, 1);
    end
    psd_block = (D2b' + D2b - 2*pass_tol*eye(p));
else
    if LMILAB
        p = dim(D2b, 1);
        Pvar = lmim('P', n, n, 'sym');
    else
        p = size(D2b, 1);
     
        Pvar = sdpvar(n, n);
    end
    
    
    psd_block = [Pvar - A2b' * Pvar * A2b,  -A2b' * Pvar * B2b + C2b';
         (-A2b' * Pvar * B2b +C2b')', (D2b' + D2b - 2*pass_tol*eye(p)) - B2b'* Pvar * B2b];
end

end

