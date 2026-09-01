%% Define LMI constraints
A=[-3 2 -1;2 -2 3;0 1 -4];
B=[0;0;1];
C=[1 0 0;2 0 1];
D=[0;1];

X=lmim('X',3,3,'sym');
ga=lmim('ga');

% observe the automatic generation of partitions
lm=[A'*X+X*A X*B C';B'*X -ga D';C D -ga*eye(2)]
see(lm)

% extract data lm.A, lm.B, lm.C lm.bl, lm.rpar, lm.cpar
[Ad,Bd,Cd,bl,rpar,cpar]=data(lm)

disp('==========================================================')

%% no automatic reduction of map representations 
see(X+X')
see(2*X)
%%
trace(lm)
trace(A'*X+X*A)+trace(-ga)+trace(-ga*eye(2))
disp('==========================================================')

%% size of array
size(lm)
% dimension of values of map 
dim(lm)
%arrays of lmim allowed
lma(1)=lm;
lma(2)=-X;
%size of array of lmim 
size(lma)
%dim(lma) not allowed for simplicity.

%most operations are not allowed for arrays.
%note that lmim(lma) leads to an error.


disp('==========================================================')

%% substitution of one variable by value
see(subs(lm,X,eye(3)))

disp('==========================================================')
%arrays of lmim are allowed in substitutions 
%to generate arrays of lmim 
subs(lma,{X,'ga'},{eye(3),3})
double(lma,{"X",ga},{eye(3),3})

%in partial substitution, non-constant maps lead to empty doubles:
subs(lma,"X",eye(3))
d=double(lma,X,eye(3))

%eig is allowed for lmival array e.g. to quickly 
% check left-hand side of LMI at solution
eig(d)

eig(lmival({rand(3,3),eye(3),rand(3,4)}))

disp('==========================================================')


