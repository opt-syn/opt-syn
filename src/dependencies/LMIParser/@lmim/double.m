function s_lmival=double(varargin)
% function s_lmival=double(varargin)
%
% Generate s_lmim=subs(varargin) and turn the value of 
% constant maps an entry (doulble) of s_lmival. For
% non-constant maps, the entry of s_lmival is empty.

p=subs(varargin{:});
s=lmival;
for j=1:numel(p);
    if numel(p(j).bl)==0
        s(j)=p(j).A;
    else
        s(j)=lmival([]);
    end
end
s_lmival=s;
end
