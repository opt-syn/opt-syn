function s=syzf(p);
%function s=syzf(p);
%
%Compute multiplier in algorithm synthesis. 
%
%The theory is exposed in C.W. Scherer, Ch. Ebenbauer, T. Holicki, 
%Optimization Algorithm Synthesis based on Integral Quadratic Constraints: A Tutorial, 
%62nd IEEE Conference on Decision and Control. 
%
%An extended version is available on arXiv under https://doi.org/10.48550/arXiv.2306.00565
%All references in the code are related to these paper.
%
%The input is collected in the structure p with following fields:
%p.l and p.L:         Strong conexity and smoothness parameters.
%p.rho:               Target convergence rate.
%p.P:                 Full plant for synthesis and simulations.
%                     Extract plant for synthesis with indices
%                     p.iz,p.iw,p.iy,p.iu.
%p.iz p.iw p.iy p.iu: Index vectors defining input-output signals of plant for synthesis. 
%p.alpha:             Characteristic polynomial of basis for multiplier (Matlab format)
%
%Results are collected in structure s, which is p with following additional fields:
%s.ga:                Computed passivity index (see below).
%                     If negative, the design is successful.
%s.X:                 Computed matrix X as in (42)
%s.Y11:               Computed matrix \tilde{Y} as in (48)
%s.mul:               Computed multiplier (34) for dimension 1 as ss object
%s.cd:                Computed value of [c d] the bold [Cf,Df] pair in paper

if ~isfield(p,'rho');p.rho=1;end;

%define plant (36) for synthesis from p.P with index vectors p.iz p.iw p.iy p.iu
Po=p.P([p.iz p.iy],[p.iw p.iu]);
Ts=Po.Ts;

%define dimension 
nz=length(p.iz);
nw=length(p.iw);
ny=length(p.iy);
nu=length(p.iu);

%define new index vectors for ease of access
iz=1:nz;
iw=1:nw;
iy=nz+(1:ny);
iu=nw+(1:nu);


if ~isfield(p, 'terminal')
    p.terminal = false;
end


%extract parameters
L=p.L;
m=p.m;
rho=p.rho;
al=p.alpha;
l=length(al)-1; %degree of polynomial = sizes of state matrix of basis 

if l>0 & al(1)~=1;error('Highest coefficient in p.al not one.');end;

%generate plant above (37) 
%rho-weighting and loop transformation 
%performed in reversed order (these operations commute)

%loop transformation of plant
I=eye(nw);
T=[-I (L-m)*I;I m*I];
Pl=lft(T,Po,nw,nw);
if norm(Pl.d(iy,iu))>0;error('Cannot allow direct feedthrough from u to y.');end;

%rho transformation plant
P=rhotrafo(Pl,rho);

%construct basis transfer matrix (see (33) and (34))
%Af   Bf
%Cfb  Dfb  (output matrices for basis) 
%
%actual filter output matrices defined as
%
%Cf=[0 c]*Cfb and Df= [d 0]*Dfb and thus (Cf Df)=(c d)
%
%with SDP variables c and d defined below

%construct realization of filter basis 
if l>0;    
    Af=jordanb(l);Bf=zeros(l,1);Bf(end)=1;
    Af(end,:)=-al(end:-1:2);
    Cfb=[zeros(1,l);eye(l)];
    Dfb=[1;zeros(l,1)];
else
    Af=[];Bf=[];Cfb=[];Dfb=1;
end;

%construct corresponding system
sfb=ss(Af,Bf,Cfb,Dfb,Ts);

%weights for multiplier set (see (32))
law=evalfr(sfb,rho);

%Filter parameters and constraints
d=sdpvar(1,1);
c=sdpvar(1,l,'full');

%Cfla and Dfla are bold Cf and Df in (33)
Cfla=[0 c]*Cfb;
Dfla=[d zeros(1,l)]*Dfb;
lmi=[];

%multiplier constraints (32)
if l>0;
    lmi=lmi+[c*ctrb(Af,Bf)<=0];
end;
lmi=lmi+[[d c]*law>=0];

%extract subsystem (Atilde, Btilde,Ctilde_z,Dtilde_z) 
%from system above (37) - change of notation to shorten formulas
[A,B,C1,E]=ssdata(P(iz,iu));
n=size(A,1);

%construct K,L,M,N (the latter two bold) in Theorem 8
if l>0;
    % K from 1)
    K=sylvester(Af,-A,-Bf*C1);
    
    % controllability matrices in 2) 
    K1=ctrb(Af,Bf);
    K2=ctrb(Af,Bf*E-K*B);
    % direct definition of inverse of L
    Li=K1*inv(K2);
    
    % compute al(Atilde) as in 3)
    alA=polyvalm(al,A);
    
    % compute second term on right hand side in 3)
    H=zeros(n);
    for nu=0:l-1;
        H=H+c(nu+1)*A^nu;
    end;
    % compute M by explicitly solving 3) with inverse of al(Atilde)
    M=d*eye(n)+H*inv(alA);
    
    %exploit linearity to compute N 
    N=zeros(n,l);
    for nu=1:l;
        N=N+c(nu)*sylvester(A,-Af,-B*Cfb(nu+1,:));        
    end;    
    %define full Tbold (as in proof) 
    %and lower block row Ttilde (as in Theorem 8)
    Tla=[Li -Li*K;N*Li  M-N*Li*K];
    Tsys=[N*Li  M-N*Li*K];
else
    Tsys = d*eye(n);
end;

%system with filter basis ((37) for filter basis)
hP=blkdiag(sfb,eye(ny))*P;
[hA,hB,hCb,hDb]=ssdata(hP);
nh=size(hA,1);

%put in parameters to generate (37)
%hA  hBw  hB
%hCz hDzw hDz
%hC  hDw     
hC=blkdiag([d c],eye(ny))*hCb;
hD=blkdiag([d c],eye(ny))*hDb;

hBw=hB(:,iw);
hCz=hC(iz,:);
hDzw=hD(iz,iw);

hC=hC(iy,:);
hDw=hD(iy,iw);

%determine annihilators in Theorem 8
U=null([hC hDw]);
nu=size(U,2);
V=null([B' E']);
nv=size(V,2);

X   = sdpvar(l+n);   % primal certificate
Y11 = sdpvar(n);     % partial dual certificate
ga  = sdpvar(1,1);   % passivity parameter

if p.terminal
    Z = sdpvar(l, l);
else
    Z = 0;
end

%bound passivity parameter from below
lmi=lmi+[ga>=-100];

%coupling (48)
lmi=lmi+[[Y11 Tsys;Tsys' X]>=0];

%matrices for primal lmi (42)
F1=[hA       hBw;
    eye(nh) zeros(nh,nw)]*U;
F2=[hCz      hDzw;
    zeros(nw,nh) eye(nw)]*U;

%replace right-lower block [0 1;1 0] in middle matrix of (42) by [0 1;1 -ga]
%ga is interpreted as a input passivity index which is minimized

%performance matrix with passivity parameter to minimize
P=[0 1;1 -ga];
lmi=lmi+[F1'*blkdiag(X,-X)*F1+F2'*P*F2<=0];

%matrices for dual lmi (48)
Fd1=V'*[-eye(n) A; zeros(nz,n) C1];
Fd2=V'*[zeros(n,nz) Tsys*hBw;
       -eye(nz) hDzw];
   
%replace right-lower block [0 1;1 0] in middle matrix of (48) by [ga 1;1 0]
%this is the inverse of [0 1;1 -ga], as motivated by dualization
   
%dual performance matrix with passivity parameter to minimize
P2=[ga 1;1 0];
lmi=lmi+[Fd1*blkdiag(Y11,-Y11)*Fd1'+Fd2*P2*Fd2'>=0];%=ga*eye(nv)

%bound Df term of multiplier 
lmi=lmi+[d<=1];

%use lmilab as solver
% opt = sdpsettings('verbose', 0, 'solver', 'lmilab');
opt = sdpsettings('verbose', 0, 'solver', 'mosek');
%minimize ga 
t = optimize(lmi, ga, opt);
s=p;
s.ga=double(ga);
s.X=double(X);
s.Y11=double(Y11);
s.mul=double([d c])*sfb;
s.cd=double([c d]);

