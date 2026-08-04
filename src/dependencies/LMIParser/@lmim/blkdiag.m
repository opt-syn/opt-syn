function s=blkdiag(varargin)
% function s=blkdiag(varargin)
%
% Block diagonal combination of lmim objects.
% Partitions inherited. 

n=length(varargin);
pin=varargin;
s=lmim(pin{1});
for j=2:n
    p=lmim(pin{j});
    
    [ks,ms]=size(s.A);
    [kp,mp]=size(p.A);
    
    if ks==0 && ms == 0
        s = p;
    elseif kp~=0 && mp ~= 0
        Z1=lmim(zeros(ks,mp));
        Z1=lmim(Z1,s.rpar,p.cpar);
        Z2=lmim(zeros(kp,ms));
        Z2=lmim(Z2,p.rpar,s.cpar);
        %          s.cpar p.cpar
        %s.rpar    s      0
        %p.rpar    0      p
        s=[s Z1;Z2 p];            
    end
end
end