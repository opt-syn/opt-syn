function so=subs(s,p,val)
% function so=subs(s,p,val)
% Evaluate maps (s.lmim) and cost (s.cost)in lmis object s to
% generate the lmim arrays lmv and cov.
%
% Input interpretation indentical to subs(s,p,val)
% for s being an lmim object.

lmv=lmim;
cov=lmim;

%if no input p and val, take variables, values from s
if nargin<2
    p=s;
end
%if no input val, take variable s, values from p
if nargin<=2
    if length(p.val)==0
        error('No values detected.')
    else
        %if values exist, evaluate p.lmim and p.cost
        if ~isempty(s.cost)
            cov=subs(s.cost,p);
        end
        for j=1:length(p.lmim)
            lmv(j)=subs(p.lmim(j),p);
        end
    end
end
if nargin==3
    %If input val exists, take variable list from p.
    %Must comply with (var,val) input of lmim/subs!
    if ~isempty(s.cost)
        cov=subs(s.cost,p,val);
    end
    for j=1:length(s.lmim)
        lmv(j)=subs(s.lmim(j),p,val);
    end
end
so=s;
so.lmim=lmv;
so.cost=cov;
so.val=[];
so.dia=[];
end
        