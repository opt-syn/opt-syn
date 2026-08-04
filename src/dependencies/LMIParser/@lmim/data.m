function [A,B,C,bl,rpar,cpar]=data(p)
% Extracts A,B,C,bl,rpar,cpar from first entry of lmim object
% array.
if numel(p)>1;
    error('Argument must not be an array.')
    %p=p(1);
end
A=p.A;B=p.B;C=p.C;
bl=p.bl;
rpar=p.rpar;
cpar=p.cpar;
end