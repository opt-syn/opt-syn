%% Define LMI constraint and LMI system 
A=[-3 2 -1;2 -2 3;0 1 -4];
B=[0;0;1];
D=[1;0];

X=lmim('X',3,3,'sym');
ga=lmim('ga',1,1);
la=lmim('la',1,2);
C=[la 3;1 0 0];

lm=[A'*X+X*A X*B C';B'*X -ga D';C D -ga*eye(2)];
lmi=lmis(lm);
lmi=lmis(lmi,1*trace(X)+ga,'c'); %adjoin cost function 
lmi1=lmisolve(lmi);
lmi1
double(X,lmi1)
double(C,lmi1)

%% to be expanded with examples