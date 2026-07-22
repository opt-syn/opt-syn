function [P] = delay_subsystem(dmax1, dmax2, r1, r2)
%DELAY_INT_SUBSYSTEM Create a linear system that performs optimization
%over delay-type communication channels.
%
%the integrator is not included.
%
%dmax1: maximum delay before the oracle
%dmax2: maximum delay after the oracle
%r1:    current delay before the oracle
%r2:    current delay after the oracle
% 
% G1 = 1/z^d1, G2 = 1/z^d2, connected with an integrator
%
%TODO: implement multi-block operator delays

Nd = dmax1+dmax2;

%store the complete state
%only tap out as appropriate
if dmax1 > 0
    % A1 = spdiags(1,1,dmax1, dmax1);
    A1 = sparse(2:dmax1, 1:dmax1-1, ones(dmax1-1, 1), dmax1, dmax1);
    B1 = zeros(Nd, 1);
    B1(1, 1) = 1;
else
    A1 = [];
    B1 = zeros(Nd, 1);
end
if dmax2 > 0
    % A2 = spdiags(1,1,dmax2, dmax2);
    A2 = sparse(2:dmax2, 1:dmax2-1, ones(dmax2-1, 1), dmax2, dmax2);
    B2 = zeros(Nd, 1);
    B2(dmax1+1, 1) = 1;
else
    A2 = [];
    B2 = zeros(Nd, 1);
end
Aint = [];

A = blkdiag(A1, A2);
B = [B2, B1];

n = size(A, 1);



% C = sparse(2, Nd, 1, 2, Nd);
C = sparse(2, Nd);
if r1 > 0   
    C(1, r1) = 1;
end

if r2 > 0   
    C(2, dmax1 + r2) = 1;
end


D = sparse(2, 2);
if r1==0
    %if no delay, directly output the gradient
    D(1, 2)=1;
end

if r2 ==0
    %if no delay, send directly to the controller
    D(2, 1)=1;
end

%assemble into a linear system
Ts = 1;
sP = sparss(A, B, C, D, 'InputName', {'w', 'u'}, 'OutputName', {'z', 'y'});
sP.Ts = Ts;

statename = cell(Nd, 1);
for i = 1:dmax1
    statename{i} = sprintf('ud%d', i);
end

for i = 1:dmax2
    statename{i+dmax1} = sprintf('wd%d', i);
end
P = full(sP);
P.StateName = statename;

end