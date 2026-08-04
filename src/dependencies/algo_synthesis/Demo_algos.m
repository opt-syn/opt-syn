%Demo algorithm synthesis based on causal OZF multipliers. 
%
%The theory is exposed in C.W. Scherer, Ch. Ebenbauer, T. Holicki, 
%Optimization Algorithm Synthesis based on Integral Quadratic Constraints: A Tutorial, 
%62nd IEEE Conference on Decision and Control. 
%
%An extended version is available on arXiv under https://doi.org/10.48550/arXiv.2306.00565
%All references in the code are related to these paper.
%
%Generation of Figs. 8-10.

clear
close all
addpath('Code')          

z=zpk('z');
%static gain 1 system
G1=ss(1);
%one step delay
Gd1=1/z;
Gd2=1/z;
%unstable system (Case 5 in Fig. 7)
Gs=ss((z-.5)/((z+1.05)*(z+.5)));

p.m=1;
p.L=10;
psim=p;

% p.alpha = [1 0 -1e-4];
p.alpha = 1;
s0=syn(Gd1,G1,p);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Case 1: Do recover gd

%static multiplier with al(z)=1 and G1(z)=G2(z) is gradient descent
p.alpha=[1];
sgr=syn(G1,G1,p)
zpk(sgr.K) %here the integrator not incluced 

%standard gradient descent with integrator 
al=2/(p.m+p.L);
A=[1];B=[-al];C=[1];
zpk(ss(A,B,C,0,-1))

%% Case 2: Do revover tm
%dynamic multiplier with alpha(z)=z (z0=0, l=1 in Step 2) of procedure in Sec. VI)
%and G1(z)=G2(z) is triple momentum 

p.alpha=[1 0];
stm=syn(G1,G1,p);
zpk(stm.K) %integrator not incluced 

%standard triple momentum with integrator
rho=1-1/sqrt(p.L/p.m);
al=(1+rho)/p.L;
be=rho^2/(2-rho);
ga=rho^2/((1+rho)*(2-rho));
A=[1+be -be;1 0];B=[-al;0];C=[1+ga -ga];
zpk(ss(A,B,C,0,-1))

%% Synthesis for Cases 3,4,5 as described in paper 
p.alpha=[1 0 -1e-4];

s1=syn(G1,Gd1,p)
s2=syn(Gd1,Gd1,p)
s3=syn(G1,Gs,p)     
save('algos')

%% Simulation 1: Standard Gradient fails for delayed plant
load ('algos')
psim.L=10;
psim.rho=s2.rho;

%use gradient descent algorithm sgr.K for plant in Case 4
psim.alg=lft(s2.P,sgr.K);

psim.col=[1 0 0]; %line color
figure(1);clf;

%simulation for quadratic function 
T=10;
ploqua(T,psim)

%simulation for non-quadratic function 
T=100;
figure(2);clf;
plonq(T,psim)

%use designed algorithm in Case 4
psim.alg=lft(s2.P,s2.K);
co=lines;psim.col=co(4,:); %line color

%simulation for quadratic function 
figure(1);
ploqua(T,psim)

%simulation for non-quadratic function 
figure(2);
plonq(T,psim)

%% Simulation 2: Case 5 with unstable plant
psim.rho=s3.rho;

%use designed algorithm from Case 5 for L=10 with function for L=11
psim.alg=lft(s3.P,s3.K);
psim.col=[1 0 0]; %line color
psim.L=11;

%simulation for quadratic function 
T=100;
figure(3);clf;
ploqua(T,psim)

%use designed algorithm from Case 5 for L=10 with function for L=10
psim.L=10;
co=lines;psim.col=co(5,:); %line color
figure(3);
ploqua(T,psim)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s=syn(G1,G2,p);

%For two systems G1 and G2 (Figure 6) and data p.m, p.L,
%set up plant for simulation, compute optimal rate and design controller.

%Generalized plant for synthesis and simulation
z=zpk('z');
Int=ss(1/(z-1));
systemnames='G1 G2 Int';
inputvar='[w;u]';
outputvar='[G1;Int]';
input_to_G1='[u]';
input_to_G2='[w]';
input_to_Int='[G2]';
p.P=sysic;

%indices of performance and control channel 
p.iw=1;p.iu=2;
p.iz=1;p.iy=2;


% K = -0.03;
% palg =lft(p.P,K);

%compute best rate by bisection over rho
s=syzfb(p);

%increase rate to guarantee success of controller synthesis 
p.rho=s.rho*1.0001;
p.mul=s.mul;

%design controller for given rate and multiplier 
s=syzfco(p);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plonq(T,p);
%% Simulation with non-quadratic function
%  p.alg algorithm as LTI system

%define function symbolically 
syms x y
fs(x,y) =p.m/2*(x^2+y^2)+(p.L-p.m)*log(exp(-x)+exp(x/3+y)+exp(x/3-y));

%compute gradient and translate to Matlab function 
fp=sym2fun(fs,p);
b=[100;-500];
p.grad=@(z) fp.g(z-b);

p.T=T;
p.d=2;
%run algsim 
s=algsim(p);

%plot results 
wnorm=sqrt(diag(s.w'*s.w));
subplot(121);
title('Iterates: Non-Quadratic')
plot(s.t,s.z,'-','Color',p.col,'MarkerFaceColor',p.col);hold on;grid on
subplot(122);
title('Norm of Gradient: Non-Quadratic')
semilogy(s.t,wnorm,'-','Color',p.col,'MarkerFaceColor',p.col);hold on;grid on

%plot rate for which synthesis has been performed
semilogy(s.t,p.rho.^(s.t),'','Color','k')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ploqua(T,p);
%% Simulation with quadratic function
U=[0 -1;1 0];
Q=U'*blkdiag(p.m,p.L)*U;
b=[100;-500];
p.grad=@(z) Q*(z-b);

p.T=T;
p.d=2;
%run algsim 
s=algsim(p);

%plot results 
wnorm=sqrt(diag(s.w'*s.w));
subplot(121);
title('Iterates: Quadratic')
plot(s.t,s.z,'-','Color',p.col,'MarkerFaceColor',p.col);hold on;grid on
subplot(122);
title('Norm of Gradient: Quadratic')
semilogy(s.t,wnorm,'-','Color',p.col,'MarkerFaceColor',p.col);hold on;grid on

%plot rate for which synthesis has been performed
semilogy(s.t,p.rho.^(s.t),'','Color','k')
end
