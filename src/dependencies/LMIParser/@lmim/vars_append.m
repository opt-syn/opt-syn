function [var_string, bl_lmibl] = vars_append(p, var, bl)
%VARS_APPEND undefined

% function [var_string,bl_lmibl]=vars_append(p_lmim)
%
% Extract unique variables from p.
%
% Generates names var and variable definitions bl of all 
% variables involved in p.

%code from lmim.vars
    
if ~isempty(p)
    for i=1:numel(p)
        %pick i-th lmim object in array
        s=p(i);

        [var, bl] = vars_append_inner(s, var, bl);
        
    end

end


var_string=var;
bl_lmibl=bl;

end