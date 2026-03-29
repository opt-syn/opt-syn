function s=ablkdiag(varargin)
% Anti-block diagonal combination of sdpss objects:
% s=ablkdiag(s1,s2,s3,...) is [0 ... 0 0 s1;0 ... s2 0;0 ... 0 s3 0 0;...]
n=length(varargin);
sinp=varargin;
s=sdpss;
for j=1:n
    s1=s;
    s2=sdpss(sinp{j});
    s.A=ablkdiag(s1.A,s2.A);
    s.B=ablkdiag(s1.B,s2.B);
    s.C=ablkdiag(s1.C,s2.C);
    s.D=ablkdiag(s1.D,s2.D);
end

function a=ablkdiag(b,c);
[n1,m1]=size(c);[n2,m2]=size(b);
a=[zeros(n2,m1) b;c zeros(n1,m2)];
end

end