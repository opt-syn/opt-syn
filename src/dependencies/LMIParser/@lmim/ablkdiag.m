function s=ablkdiag(varargin)
% function s=ablkdiag(varargin)
%
% Anti-block-diagonal combination of lmim objects.
% Starting from upper right to lower left.
% Partitions inherited. 


n=length(varargin);
pin=varargin;
s=lmim(pin{1});
for j=2:n
    p=lmim(pin{j});
    
    [ks,ms]=size(s.A);
    [kp,mp]=size(p.A);
    

    Z1=lmim(zeros(ks,mp));
    Z1=lmim(Z1,s.rpar,p.cpar);
    Z2=lmim(zeros(kp,ms));
    Z2=lmim(Z2,p.rpar,s.cpar);

    %          mp:p.cpar ms:s.cpar
    %ks:s.rpar    0      s
    %kp:p.rpar    p      0
    s=[Z1 s;   
       p Z2];
end
end