function [sys_lti] = periodic_lift(sys)
%PERIODIC_LIFT lift a periodic linear system into an LTI system



Nss = length(sys);

A = cell(1, Nss);
B = cell(1, Nss);
C = cell(1, Nss);
D = cell(1, Nss);


% if any(cellfun(@isempty, sys))
%     sys_lti = [];
% else
    
    for i = 1:Nss
        [A{i}, B{i}, C{i}, D{i}] = ssdata(sys{i});
    end
    
    [n, m] = size(B{1});
    p = size(C{1}, 2);
    
    
    Aall = eye(n);
    for i = 1:Nss
        Aall = A{i}*Aall;
    end
    
    Ball = B{end};
    AfundB = eye(n);
    for i = (Nss-1):-1:1
        AfundB = AfundB * A{i+1};
        Ball = [AfundB * B{i}, Ball];
    end
    
    Call = C{1};
    AfundC = eye(n);
    for i = 2:Nss
        AfundC = A{i-1} * AfundC;
        Call = [Call; C{i} * AfundC];
    end
    
    Dall = blkdiag(D{:});
    
    for i = 2:Nss
        for j = 1:i-1
    
            ind_i = (1:m) + (i-1)*m;
            ind_j = (1:m) + (j-1)*m;
    
            AfundBC = eye(n);
            for k = (j+1):(i-1)
                AfundBC  = A{k} * AfundBC;
            end
            Dall(ind_i, ind_j) = C{i} * AfundBC * B{j};
        end
    end
    
    sys_lti = ss(Aall, Ball, Call, Dall, 1);
    
    end
% end