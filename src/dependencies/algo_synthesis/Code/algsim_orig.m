function s=algsim_orig(p);
%function s=algsim(p);
%
%Simulate closed loop system.
%
%The theory is exposed in C.W. Scherer, Ch. Ebenbauer, T. Holicki, 
%Optimization Algorithm Synthesis based on Integral Quadratic Constraints: A Tutorial, 
%62nd IEEE Conference on Decision and Control. 
%
%An extended version is available on arXiv under https://doi.org/10.48550/arXiv.2306.00565
%All references in the code are related to these paper.
%
%Generate trajectories through recursion (Figure 1)
%
%x_(t+1)=Ax_t+Bw_t, w_t=nabla f(z_t)
%z_t    =Cx_t+Dw_t 
%
%with initial condition zero. A nonzero trajector is 
%generated with a function f for which nabla f(0)<>0.
%
%The input is collected in the structure p with following fields:
%p.alg            Closed loop system in ss format as used for simulation 
%p.d              Dimension for Kroneckering 
%p.grad           Matlab function to compute gradient 
%
%%Results are collected in structure s with following additional fields:
%s.x              state trajectory
%s.w              input trajectory
%s.z              output trajectory
%s.t              time instants

[a,b,c,d]=ssdata(p.alg);
if norm(d(1,1))>1e-6;
    error('Direct feedthrough not zero.');
end;

A=kron(a,eye(p.d));
B=kron(b(:,1),eye(p.d));
D=kron(d(:,1),eye(p.d));
C=kron(c,eye(p.d));

n=size(A,1);

if isfield(p,'x');
    x=kron(p.x,ones(p.d,1));
else
    x=zeros(n,1);
end
z=C(1:p.d,:)*x;
w=p.grad(z);
xv=[x];
wv=[w];
zv=[C*x+D*w];
for t=1:p.T-1;
    w=p.grad(C(1:p.d,:)*xv(:,end));
    xp=A*xv(:,end)+B*w;
    zp=C*xv(:,end)+D*w;
    wp=p.grad(zp(1:p.d));
    xv=[xv xp];
    zv=[zv zp];
    wv=[wv wp];
end;
s.x=xv;
s.w=wv;
s.z=zv;
s.t=0:p.T-1;
