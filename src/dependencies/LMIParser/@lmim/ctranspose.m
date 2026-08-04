function s=ctranspose(p)
% function s=ctranspose(p)
%
% Transpose lmim object p. 
% From A+B*X*C to A'+B'*X'*C'.

p=lmim(p);
s=p;
s.A=p.A';
s.B=p.C';
s.C=p.B';
s.rpar=p.cpar;
s.cpar=p.rpar;
h=p.bl;
for j=1:length(h);
    %use transposition operator for lmiblk
    h(j)=h(j)';
end
s.bl=h;
% for j=1:length(p.bl)
%     %use transposition operator for lmiblk
%     s.bl(j)=p.bl(j)';
% end
