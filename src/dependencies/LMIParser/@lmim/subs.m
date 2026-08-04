function s_lmim=subs(p_lmim,var,val)
% function so_lmim=subs(p_lmim,var,val)
% Evalute maps p_lmim in variables defined by var and val
%
% If var is lmis: Ignore val, values taken from fields (p.val,p.var)
%
% If var is list of variables (see function exvar) then
% val is corresponding list of values (see function checkval).
%
% Ordering of var must match ordering of val.

s=p_lmim;
if nargin==1;
    so=p_lmim;    
else
    if nargin==2
        %val ignored!        
        p=lmis(var);
        var=p.var;
        val=p.val;
    end
    so=lmim;
    for j=1:numel(s);
        so(j)=subsloc(s(j),var,val);
    end
end
s_lmim=so;
end
