function s=syzfco(p);
%function s=synco(p);
%
%Compute controller/algorithm for given multiplier in algorithm synthesis. 
%
%The theory is exposed in C.W. Scherer, Ch. Ebenbauer, T. Holicki, 
%Optimization Algorithm Synthesis based on Integral Quadratic Constraints: A Tutorial, 
%62nd IEEE Conference on Decision and Control. 
%
%An extended version is available on arXiv under https://doi.org/10.48550/arXiv.2306.00565
%All references in the code are related to these paper.
%
%For given multiplier, transform problem to H-infty design 
%and construct controller with existing hinfsyn routine in Matlab.
%
%p.l and p.L:         Strong conexity and smoothness parameters
%p.rho:               Target convergence rate
%p.P:                 Complete plant. Extract plant for synthesis with indices.
%p.iz p.iw p.iy p.iu: Indizes of input-output signals of plant.
%p.mul:               Multiplier in ss format e.g. as computed by syzf.
%
%Results are collected in structure s, which is p with following additional fields:
%
%s.order_plant_full_red: Array with orders of [plant, controller, reduced controller]
%s.pa_ga_full:           Array with [input passivity index of closed loop,H-infinty norm of transformed closed loop]
%s.pa_ga_red:            Same for reduced controller.
%                        These numbers must be [<0 and <1] for success.
%s.Kfull:                Full controller according to (40)
%s.K:                    Reduced controller according to (40)

%extract plant and keep sample time
Po=p.P([p.iz p.iy],[p.iw p.iu]);
Ts=Po.Ts;

%the next parts are copied from the beginning of syzf
nz=length(p.iz);
nw=length(p.iw);
ny=length(p.iy);
nu=length(p.iu);

iz=1:nz;
iw=1:nw;
iy=nz+(1:ny);
iu=nw+(1:nu);

L=p.L;
m=p.m;
rho=p.rho;

I=eye(nw);
T=[-I (L-m)*I;I m*I];
Pl=lft(T,Po,nw,nw);

P=rhotrafo(Pl,rho);

%construct plant with multiplier for negative real synthesis
sys=blkdiag(p.mul,eye(ny))*P;

%transform plant to design controller by H-infinty synthesis
%involves loop transformation to change negative realness to gain bound one
E=eye(nz);
T=[E sqrt(2)*E;sqrt(2)*E E];
syshi=lft(T,sys,nz,nz);

%design controller with standard routine in Matlab
opt=hinfsynOptions;
opt.Method='lmi';
opt.RelTol=1e-7;
[Kfull,cl,ga,info] = hinfsyn(syshi,ny,nu,opt);

%reduce controller by truncation of modes
K=minreal(Kfull,1e-5);

%check for full and reduced controller whether specs are achieved
%do that for both negative realness and gain bound one
s=p;
%return order of plant, full controller, reduced controller
s.order_plant_full_red=[size(syshi.a,1),size(Kfull.a,1),size(K.a,1)];

%return passivity indices and gain to verify success 
s.pa_ga_full=[-getPassiveIndex(-lft(sys,Kfull),'in') norm(lft(syshi,Kfull),inf)];
s.pa_ga_red=[-getPassiveIndex(-lft(sys,K),'in') norm(lft(syshi,K),inf)];

if s.pa_ga_red(2)>=1
    % warning('Full controller design not succesful!');
    s.K = [];
else
    s.K=rhotrafo(K,1/rho);
end

if s.pa_ga_full(2)>=1
    % warning('Reduced controller design not succesful!');
    s.Kfull = [];
else
    s.Kfull=rhotrafo(Kfull,1/rho);
end

%full and reduced controller transformed with rho as in (40)

