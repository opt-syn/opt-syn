%A simple example to compute the H-infinity norm of a system:

A=[-3 2 -1;2 -3 3;0 1 -4];
B=[0;0;1];
C=[1 0 0;0 1 0];
D=[0;1];
n=size(A,1);

%define KYP variable X (the map 0+I*X*I)
X=lmim('X',n,n,'sym');

%define variable ga (once repeated)
ga=lmim('ga',1,1);

%define LMI (by default to be read as lm<=0)
lm=[A'*X+X*A X*B C';B'*X -ga D';C D -ga*eye(2)]
dim(lm)

%define cost
cost=ga;

%collect first LMI l1 into LMI system lmi
lmi=lmis(lm);

%add the cost
lmi=lmis(lmi,cost,'c');
see(lmi)

%solve the LMI with LMIlab (second argument can be LMIlab options opt)
[lmi1,info1]=lmisolve(lmi);
lmi1

%evaluate lmi constraints and cost
[lv1,cv1]=double(lmi1)

%solve the LMI via Yalmip (second argmument can be sdpsettings options)
[lmi2,info2]=yalsolve(lmi);
lmi2



%%
doc lmibl
doc lmim
doc lmis
doc lmival

