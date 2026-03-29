%% Testing of variable construction and allowed concatenations

%define scalar variable and see difference as variable and map
xv=lmibl('x') 
x=lmim('x')
see(x)
%there is no need to work with lmibl!

xr=x*eye(3)   %allowed, type changed
%xr=xv*eye(3)  %not allowed
disp('==========================================================')
%%
%Concatenation of mpas only posible if dimensions fit and
%the involved variables are compatible if they have identical names.

X=lmim('X',2,4,'full') 
Y=lmim('Y',2,2,'sym')
X2=lmim('X',2,2,'sym')

[X X]  %is fine
[X Y]  %is fine
%[X X2] %not ok
disp('==========================================================')

%% Variable construction and pitfalls
% Apart from horizontal concatenation (for combining variables),
% no further operations for lmibl objects are allowed.
% Rather perform operations for the corresponding lmim objects!
X=lmibl('X',3,3,'sym')
Y=lmibl('Y',3,3,'full');

% this is an lmibl array
Z=[X Y]
% [X;Y] %not allowed
% lmim(Z) %not allowed

Xm=lmim(X);
Ym=lmim(Y);
[Xm Ym]
[var,bl]=vars([Xm;Ym]);
var
bl(1)
bl(2)

Z=lmim(lmibl('X',3,3,'full'))
% Z=[Xm Z] %not allowed

%Partitioning
lmim([Xm;Ym],[2 3 1],[1 1 1])
disp('==========================================================')
%%
%different ways of diagonal augmentation of scalar variable
ga=lmim('ga')
Ga=drep(ga,3)
lmim('ga',3,0)
lmim(ga*eye(3))
disp('==========================================================')
%%
X=lmim('X',2,2,'full')
Y=lmim('Y',3,3,'sym');
%define with row and column partition 
Y=lmim(Y,[1 2],[1 2])
disp('==========================================================')
%% check partitions
blkdiag(X,Y)
ablkdiag(X,Y)
z=trace(Y);
see(z)
Z=blkdiag(X,z)
disp('==========================================================')
%% lmim arrays are allowed
Y(2)=lmim('X',3,3)
[Y(1);Y(2)]
[Y(1) Y(2)]

%sel selects subblocks
sel(Y(2),[1 2],[1 2])
%Y([1,2],[1 2]) does not work
disp('==========================================================')

%% affine map
% constant maps
A=[1;3;1];
B=lmim([1 2 -1;2 1 3;0 1 0]);
C=lmim([0;0;1])

% affine  maps 
D=A+B*Z*C
see(A+B*Ga*C)
see(A+B*ga*eye(3)*C)
disp('==========================================================')
%% Additon and multiplication 
E=A*C'
E+Z
-Z
E-Z
%E+(-1)*Z does not work: scalar times lmim not defined

X=lmim('X',1,2);
L1=[1 2]+2*X*[1 0;3 0]
L2=[1;2]+[0;1]*X*[1;2]
see(L1*L2)
%X*[1 0;3 1]*[0;1]*X not allowed
disp('==========================================================')

%% Extract variable names 
x=lmim('x',2,1);
ga=lmim('ga');
de=lmibl('de',1,2);

exvar({de,ga*eye(3),[x;x],'X',"Y"})
[v,bl]=vars([x;ga])
bl(1)
bl(2)
%exvar([x;ga]) %not allowed

disp('==========================================================')
%% List of combinations 
% fine
blkdiag(lmim('X',1,1,'rep'),lmim('X',1,1,'rep'))
blkdiag(lmim('X',1,1,'rep'),lmim('X',2,2,'rep'))
blkdiag(lmim('X',1,1,'rep'),lmim('X',1,1,'full'))
blkdiag(lmim('X',1,1,'rep'),lmim('X',1,1,'sym'))

blkdiag(lmim('X',2,2,'rep'),lmim('X',1,1,'rep'))
blkdiag(lmim('X',2,2,'rep'),lmim('X',2,2,'rep'))
blkdiag(lmim('X',2,2,'rep'),lmim('X',1,1,'full'))
blkdiag(lmim('X',2,2,'rep'),lmim('X',1,1,'sym'))

blkdiag(lmim('X',1,1,'sym'),lmim('X',1,1,'rep'))
blkdiag(lmim('X',1,1,'sym'),lmim('X',2,2,'rep'))
blkdiag(lmim('X',1,1,'sym'),lmim('X',1,1,'full'))
blkdiag(lmim('X',1,1,'sym'),lmim('X',1,1,'sym'))

blkdiag(lmim('X',1,1,'full'),lmim('X',1,1,'rep'))
blkdiag(lmim('X',1,1,'full'),lmim('X',2,2,'rep'))
blkdiag(lmim('X',1,1,'full'),lmim('X',1,1,'full'))
blkdiag(lmim('X',1,1,'full'),lmim('X',1,1,'sym'))

blkdiag(lmim('X',2,2,'sym'),lmim('X',2,2,'sym'))
blkdiag(lmim('X',2,2,'full'),lmim('X',2,2,'full'))
blkdiag(lmim('X',2,3,'full'),lmim('X',2,3,'full')')

%% Not fine
blkdiag(lmim('X',2,2,'sym'),lmim('X',1,1,'rep'))
blkdiag(lmim('X',2,2,'full'),lmim('X',1,1,'rep'))

blkdiag(lmim('X',2,2,'sym'),lmim('X',2,2,'rep'))
blkdiag(lmim('X',2,2,'full'),lmim('X',2,2,'rep'))

blkdiag(lmim('X',2,2,'sym'),lmim('X',2,2,'full'))

blkdiag(lmim('X',2,2,'sym'),lmim('X',3,3,'sym'))
blkdiag(lmim('X',2,2,'sym'),lmim('X',2,2,'full'))

%% Not fine
lmim('X',2,3,'full')
lmim(lmibl('X',3,2,'full')')
blkdiag(lmim('X',2,3,'full'),lmim(lmibl('X',3,2,'full')'))
