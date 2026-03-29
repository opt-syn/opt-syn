z=zpk('z');
sy=ss([1/(z-1) 1/z;1/(z^2+1) z/(z^3-2)]);
G=sdpss(sy)
F=sdpss(sdpvar(2,2))
H=sdpss(rand(2,2))
%operations fine if first is sdpss object 
G*F
%does not work otherwisem
% F+H+G

%% check transfer function evaluation in matrix
M=rand(3,3)-rand(3,3);
%M=[1 2 -2;3 1 -2;1 1 2];
%M=diag([3,4]);
z=zpk('z')
al=(z-1)^2*(z+2);
be=(z+3)*(z+1);

pal=tf(al).Numerator{1};
pbe=tf(be).Numerator{1};

R1=polyvalm(pbe,M)*inv(polyvalm(pal,M))
sy=sdpss(ss(be/al));

R2=evalm(M,sy)
R1-R2

%% check sylvester solver
        M=rand(3,2)-rand(3,2);
A=rand(3,3)-rand(3,3);
B=rand(2,2)-rand(2,2);
X=sylvester(A,B,M);
Y=syl(A,B,M);
X-Y

%%

G=sdpss(rss(3,2,2));
E=sdpss(sdpvar(2,2));
H=sdpss(rand(4,2));

sy=[H;G;E;E]
%%
G=sdpss(rss(3,2,2));
H=rss(3,2,2);
s=G+H
s.sys


%%
F=sdpss(1)
G=sdpss(1)
H=sdpss(1)

%%
s1=[G,E];
s2=ssostackv(s1,H)

ssovar(s1)
ssovar(s2)

mysee(s2.sys)
ssosee(s2);
%%
[s,sc,sv,v]=ssosee(s2)
s
%ssostackv(G,H)
%%
A=rand(3,3);
B=rand(3,1);
C=sdpvar(1,3);
D=1;
sy=sso(A,B,C,D)
ssosee(sy)
