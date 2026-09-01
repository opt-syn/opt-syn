function [lmv,cov]=double(varargin)
%as for subs, but returns mival array lmv, cov with value of LMI
%constraints, cost

if isnumeric(varargin{1})
    lmv = varargin{1};
    cov = [];
else
    s=subs(varargin{:});
    lmv=double(s.lmim);
    cov=double(s.cost);
end
end
        