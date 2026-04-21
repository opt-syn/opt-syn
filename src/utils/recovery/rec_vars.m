function [vrec] = rec_vars(vars, lmi_out)
%REC_VARS recover variables from a problem
%   1 argument: yalmip (use value)
%   2 arguments: lmilab parser (use double)

if nargin < 2 || (isnumeric(lmi_out) && lmi_out == 0)
    LMILAB = false;
    lmi_out = 0;
else
    LMILAB = true;
end

vrec = vars;
field = fieldnames(vars);
for i = 1:length(field)
    fc = getfield(vars, field{i});


    if isstruct(fc)
        %recursively go through structures
        fcrec = rec_vars(fc, lmi_out);
        vrec = setfield(vrec, field{i}, fcrec);
    else
        if ~isnumeric(fc)
            if iscell(fc)
                if LMILAB
                    fcrec = cell(size(fc));
                    for j = 1:length(fc)
                        if isnumeric(fc{j})
                            fcrec{j} = fc{j};
                        elseif isstruct(fc{j})
                            fcrec{j} = rec_vars(fc{j}, lmi_out);
                        else
                            fcrec{j} = double(double(fc{j}, lmi_out));
                        end
                    end
                    % fcrec = cellfun(@(c) double(double(c, lmi_out)), fc, 'UniformOutput', false);
                else
                    fcrec = cellfun(@value, fc, 'UniformOutput', false);
                end
            else
                % vrec = setfield(vrec, field{i}, value(fc));
                % fcrec = value(fc);
                if LMILAB
                    fcrec = double(double(fc, lmi_out));
                else
                    fcrec = value(fc);
                end
            end
            vrec = setfield(vrec, field{i}, fcrec);
        end
    end

    
end

end

