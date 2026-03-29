function [var_string,bl_lmibl]=vars(p_lmim)
% function [var_string,bl_lmibl]=vars(p_lmim)
%
% Extract unique variables from p.
%
% Generates names var and variable definitions bl of all 
% variables involved in p.

p=p_lmim;
var=string([]);
bl=lmibl([]);
for i=1:numel(p)
    %pick i-th lmim object in array
    s=p(i);
    for j=1:numel(s.bl)
        va=s.bl(j);
        %test whether variable va in var
        %if yes then iv is index of current variable in var
        [ind,iv]=ismember(va.na,var);

        %only add variables which are NOT transposed
        %if variable is transposed, select non-transposed version

        %if va.na not in list var
        if ind==0
            %add na to list var
            var=[var string(va.na)];
            %new length of list
            le=length(var);
            %add original (not transposed) variable to so.varbl
            if va.tr
                bl(le)=va';
            else
                bl(le)=va;
            end
        end
    end
end
var_string=var;
bl_lmibl=bl;
end
