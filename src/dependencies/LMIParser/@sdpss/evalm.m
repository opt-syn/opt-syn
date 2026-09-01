function N=evalm(M,s)
%evaluation of sdpss system s in double matrix M
%resulting in sdpvar matrix N
%
%C.W. Scherer

s=sdpss(s);
if ~isa(M,'double')
    error('First entry must be a matrix.')
end
if size(M,1)~=size(M,2)
    error('First entry must be square.')
end
if isa(s.A,'sdpvar')
    error('Property A of second entry cannot be an sdpvar.')
end

m=size(M,1);
I=eye(m);
n=size(s.A);
Ae=kron(M,eye(n))-kron(eye(m),s.A);
if min(svd(Ae))<eps
    error('Well-posedness of LFT not guaranteed.')
end
N=kron(I,s.D)+kron(I,s.C)*inv(Ae)*kron(I,s.B);
end
