function [var, bl] = vars_append_inner(s, var, bl)
%VARS_APPEND_INNER undefined
%   undefined
    for j=1:numel(s.bl)
        na=s.bl(j).na;
        %test whether variable va in var
        %if yes then iv is index of current variable in var
        if isempty(var)
            ind = 0;
            % iv = 0;
        else
            ind = any(strcmp(na,var));
            % [ind,iv]=ismember(va.na,var);
        end
    
        %only add variables which are NOT transposed
        %if variable is transposed, select non-transposed version
    
        %if va.na not in list var
        if ind==0
            %add na to list var
            var=[var string(na)];
            %new length of list
            le=length(var);

            
            %add original (not transposed) variable to so.varbl

            if isnumeric(bl) && isempty(bl)
                bl=lmibl([]);
            end
            va = s.bl(j);
            if va.tr
                bl(le)=va';
            else
                bl(le)=va;
            end
        end
    end
end