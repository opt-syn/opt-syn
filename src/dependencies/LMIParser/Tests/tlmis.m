%% Define LMI constraints
A=[-3 2 -1;2 -2 3;0 1 -4];
B=[0;0;1];
C=[1 0 0;2 0 1];
D=[0;1];

X=lmim('X',3,3,'sym');
ga=lmim('ga');

% observe the automatic generation of partitions
lm=[A'*X+X*A X*B C';B'*X -ga D';C D -ga*eye(2)]

%% Define lmi system and cost
%generate LMI system with one LMI lm<=0
lmi=lmis(lm);

%add cost (scalar-valued lmim object) 
lmi=lmis(lmi,ga,'c')

%% Solve with LMIlab and evaluate 
% solve with LMIlab
[lmi1,info1]=lmisolve(lmi,lmisolveopt)

%subsitute values, all values in lmi1 are deleted.
lmi2=subs(lmi1)

%generate lmival of values (for constant maps) after substituion 
[l1,c1]=double(lmi1) 
eig(l1)
%double of lmival array generates cell array of doubles 
max(double(eig(l1)))
lmi1.dia(2)
%% Solve with Yalmip and evaluate 
% solve with Yalmip
[lmi2,info2]=yalsolve(lmi)
l2=double(double(lmi2))
double(l1)

%% Partial substitution and solve new LMI problem 
% fix ga to value and minimize trace (smallest solution)
lma=subs(lmi,ga,12);
lma=lmis(lma,trace(X),'c');
lma1=lmisolve(lma)
lma2=yalsolve(lma)
double(double(lma1))-double(double(lma2))

