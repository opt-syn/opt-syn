function s=uminus(p)
% function s=uminus(p)
%
% Negative of lmim object (Multiplication by -1)
p=lmim(p);
s=p;
s.A=-p.A;
s.B=-p.B;
end
