function [A,B,C]=coelist(p)
% Genrates matrix coefficient lists A{j,k} and B{i,j}, C{j,k}
% j refers to index of variable block
% i,k refer to index of partition in p.rpar, p.cpar
%
% Subplock {i,k} in partition of A+B1*X1*C1 + B2*X2*C2 + ... given by
% A{i,k} + B{i,1}*X1*C{1,k} + B{i,2}*X2*C{2,k}+ ...

p=lmim(p);

jr0=0;
jc0=0;
for j=1:length(p.bl)
    jrd=p.bl(j).di(1);
    jcd=p.bl(j).di(2);
    i0=0;
    for i=1:length(p.rpar)
        id=p.rpar(i);
        B{i,j}=p.B(i0+(1:id),jr0+(1:jrd));
        i0=i0+id;
    end
    k0=0;
    for k=1:length(p.cpar)
        kd=p.cpar(k);
        C{j,k}=p.C(jc0+(1:jcd),k0+(1:kd));
        k0=k0+kd;
    end
    jr0=jr0+jrd;
    jc0=jc0+jcd;
end
i0=0;
for i=1:length(p.rpar)
    k0=0;
    for k=1:length(p.cpar)
        id=p.rpar(i);
        kd=p.cpar(k);
        A{i,k}=p.A(i0+(1:id),k0+(1:kd));
        k0=k0+kd;
    end
    i0=i0+id;
end
end
