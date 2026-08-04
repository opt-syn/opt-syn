function var=exvar(p)
% function var=extvar(p)
% Extracts string array var of names of all variables in p.
%
% p  : string, char, lmim, lmibl object
%    : lmim and lmibl arrays can involve only ONE variable!
% p  : cell of string, char, lmim, lmibl objects (mixtures allowed!)
%    : lmim and lmibl entries (possibly arrays) can invovle only ONE variable!
% 
% Serves to easily generate list of variable names for substituion.
% Such lists can as well be composed with the command vars. 

%no output if p is empty
if isempty(p)
    var=[];
    return
end

if ~iscell(p)
    p={p};
end

st=[];
for j=1:numel(p)
    s=p{j};
    switch class(s)
        case 'char'
            st=[st string(s)];
        case 'string'
            st=[st string(s)];
        case {'lmibl','lmim'}
            h=vars(s);
            if numel(h)>1
                error(['An lmim element comprises more than one variable.'])
            else
                st=[st h];
            end
        otherwise
            error('Variable list contains objects not allowed.')
    end
end
%if unique(st)
if numel(unique(st))<numel(st)
    st
    error('Generated list of variable names has identical entries.')
end
var=st;
