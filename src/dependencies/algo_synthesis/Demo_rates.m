%Demo algorithm synthesis based on causal OZF multipliers.
%
%The theory is exposed in C.W. Scherer, Ch. Ebenbauer, T. Holicki, 
%Optimization Algorithm Synthesis based on Integral Quadratic Constraints: A Tutorial, 
%62nd IEEE Conference on Decision and Control. 
%
%An extended version is available on arXiv under https://doi.org/10.48550/arXiv.2306.00565
%All references in the code are related to these paper.
%
%Generation of Fig. 7 with optimal rates for synthesis.

clear
addpath('Code')
z=zpk('z');

%static gain 1 system
G1=ss(1);
%one step delay
Gd=1/z;
%unstable system (Case 5 in Fig. 7)
Gs=ss((z-.5)/((z+1.05)*(z+.5)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute rates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%grid for values of L with m=1
p.Lv=logspace(log10(10),log10(200),20);

%static multiplier and G1(z)=G2(z) is gradient descent
p.alpha=[1];
rhovgd=rates(G1,G1,p)

%dynamic multiplier with alpha(z)=z (z0=0, l=1 in Step 2) of procedure in Sec. VI)
%and G1(z)=G2(z) is triple momentum 
p.alpha=[1 0];
rhovtm=rates(G1,G1,p)

%dynamic multiplier with alpha(z)=z^2-z0^2 (z0=0.01, l=2 in Step 2) of procedure in Sec. VI)
%generates curves of Case 3, 4, 5
p.alpha=[1 0 -1e-4];
rhov3=rates(G1,Gd,p)
rhov4=rates(Gd,Gd,p)
rhov5=rates(G1,Gs,p)

save('rates')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot rates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clf
figure(WindowState="maximized")
load('rates')
rhovg=(p.Lv-1)./(p.Lv+1);
rhovt=1-1./sqrt(p.Lv);
co=lines;
col=co(1,:);
semilogx(p.Lv,rhovg,'-o','LineWidth',1,'Color',col,'MarkerFaceColor',col);hold on;grid on
col=co(2,:);
semilogx(p.Lv,rhovt,'-o','LineWidth',1,'Color',col,'MarkerFaceColor',col);hold on;grid on
col=co(3,:);
semilogx(p.Lv,rhov3,'-o','LineWidth',1,'Color',col,'MarkerFaceColor',col);hold on;grid on
col=co(4,:);
semilogx(p.Lv,rhov4,'-o','LineWidth',1,'Color',col,'MarkerFaceColor',col);hold on;grid on
col=co(5,:);
semilogx(p.Lv,rhov5,'-o','LineWidth',1,'Color',col,'MarkerFaceColor',col);hold on;grid on
xlabel('L')
ylabel('Optimal \rho_*')
legend('\rm Case 1: Gradient decent','\rm Case 2: Triple momentum algorithm','Case 3: $G_1(z)=1$, $G_2(z)=\frac{1}{z}$','Case 4: $G_1(z)=\frac{1}{z}$, $G_2(z)=\frac{1}{z}$','Case 5: $G_1(z)=1$, $G_2(z)=\frac{z-0.5}{(z+0.5)(z+1.05)}$','Location','se','interpreter','latex','Fontsize',25)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rhov=rates(G1,G2,p);

%For two systems G1 and G2 (Figure 6) and data p.Lv,
%set up plant for synthesis and compute optimal rates for
%m=1 and for L in the list of values in p.Lv

%Generalized plant for synthesis (and simulation)
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

%strong convexity parameter
p.m=1;

%compute optimal rates for values of L in p.Lv
rhov=[];
for L=p.Lv;
    p.L=L;
    s=syzfb(p);
    rhov=[rhov s.rho]
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

