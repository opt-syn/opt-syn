function so=drep(s,d)
% function so=drep(s,d)
%
% d times repetition of s corresponds to 
% blkdiag(s,...,s) d-times

s=lmim(s);
if nargin>1 & d>1
    so=s;
    for j=2:d
        so=blkdiag(so,s);
    end
end
end
