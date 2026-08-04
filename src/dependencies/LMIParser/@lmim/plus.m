function s=plus(s1,s2)
% function s=plus(s1,s2)
%
% Sum of two lmim objects s1 and s2
% (A1+B1*X1*C1) + (A2+B2*X2*C2)=
% (A1+A2) + [B1 B2][X1 0;0 X2][C1;C2]=
% (A1+A2) + [B1*X1 B2*X2][C1;C2]

s1=lmim(s1);
s2=lmim(s2);
if ~isequal(dim(s1),dim(s2))
    error('Dimensions not compatible.')
end

s=lmim;
s.A=s1.A+s2.A;
s.B=[s1.B s2.B];
s.C=[s1.C;s2.C];


% combine row and column partitions
s.rpar=parcomb(s1.rpar,s2.rpar);
s.cpar=parcomb(s1.cpar,s2.cpar);

%put at end for validation 
% concatenation of two block cells is diagonal combination 
s.bl=[s1.bl s2.bl];

