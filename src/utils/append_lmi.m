function [lmi_new] = append_lmi(lmi_old,lmi_add, LMILAB)
%APPEND_LMI add a new LMI to a list of LMIs

%default to lmi_add >= 0
%YALMIP: this is ok
%LMIPARSER: need to flip the sign to produce -lmi_add <= 0.

if isnumeric(lmi_add)   
    lmi_new = lmi_old;
else
    if LMILAB
        if isempty(lmi_old)
            lmi_new = lmis(-lmi_add);
        else
            lmi_new = lmis(lmi_old, -lmi_add);
        end   
    else
        lmi_new = [lmi_old; lmi_add >= 0];
    end
end

end

