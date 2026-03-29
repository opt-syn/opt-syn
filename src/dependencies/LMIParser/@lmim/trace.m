function s=trace(p)
% function s=trace(p)
%
% Trace of lmim object p. 
% No partition of output. 

p=lmim(p);
[k,m]=size(p.A);
if k~=m
    error('Map does not have square values.')
else
    s=sel(p,1,1);
        for j=2:k;
            s=s+sel(p,j,j);
        end
end
end
