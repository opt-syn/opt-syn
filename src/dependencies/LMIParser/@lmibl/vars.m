function var_string=vars(p_lmibl)
% function var_string=vars(p_lmibl)
%
% Generates var of the names of 
% all variables involved in p.

p=p_lmibl;
var=[];
bl=[];
lm=lmim;
for i=1:numel(p)
    %pick i-th lmibl object in array
    lm(i)=p(i);
end
var_string=vars(lm);
end
