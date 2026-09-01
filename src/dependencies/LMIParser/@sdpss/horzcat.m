function s=horzcat(varargin)
% Horizontal concatenation of spdss objects.

n=length(varargin);
sinp=varargin;
s=sdpss;
for j=1:n
    s1=s;
    s2=sdpss(sinp{j});
    s.A=blkdiag(s1.A,s2.A);
    s.B=blkdiag(s1.B,s2.B);
    s.C=[s1.C s2.C];
    s.D=[s1.D s2.D];
end
end