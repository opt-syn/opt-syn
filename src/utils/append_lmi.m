function [lmi_new] = append_lmi(lmi_old,lmi_add, LMILAB)
%APPEND_LMI add a new LMI to a list of LMIs

%default to lmi_add >= 0
%YALMIP: this is ok
%LMIPARSER: need to flip the sign to produce -lmi_add <= 0.

if isnumeric(lmi_add)   
    lmi_new = lmi_old;
else
    if LMILAB

        if iscell(lmi_add)
                % lmi_neg = lmi_add;
                
            lmi_neg = cellfun(@(c) -c, lmi_add, 'UniformOutput', false);
        else
            lmi_neg = -lmi_add;
        end


        if isempty(lmi_old)
            lmi_new = lmis(lmi_neg);
        else
            
            lmi_new = lmis(lmi_old, lmi_neg);
        end   
    else
        lmi_new = [lmi_old; lmi_add >= 0];
    end
end

end

