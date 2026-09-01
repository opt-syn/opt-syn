function s=vertcat(varargin)
% Vertical concatenation of spdss objects.

n=length(varargin);
sinp=varargin;
s=sdpss;
for j=1:n
    s1=s;
    s2=sdpss(sinp{j});

    [ks,ms]=size(s1.D);
    [kp,mp]=size(s2.D);
    if ks==0 || ms ==0
        s = s2;
    elseif kp==0 || mp ==0
        s = s1;
    else
        s.A=blkdiag(s1.A,s2.A);
        s.B=[s1.B s2.B];
        s.C=blkdiag(s1.C,s2.C);    
        s.D=[s1.D;s2.D];
    end
end
end