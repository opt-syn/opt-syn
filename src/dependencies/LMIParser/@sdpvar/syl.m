function X = syl(A, B, C)
%SYLVESTER Solve Sylvester Equation.
%   X = SYL(A,B,C) solves the Sylvester equation A*X + X*B = C,
%   where A is a m-by-m matrix, B is a n-by-n matrix, and C
%   is an m-by-n matrix or sdpvar.
%   
%   X is an m-by-n matrix or sdpvar.
%   The implementation might not be efficient.
%
%C.W. Scherer

n=size(A,1);
m=size(B,1);
Ae=kron(eye(m),A)+kron(B',eye(n));
Xe=Ae\reshape(C,n*m,1);
X=reshape(Xe,n,m);