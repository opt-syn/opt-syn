function  [di,nc]=dim(p,dim)
% function  [di,nc]=dim(p,dim)
% Displays or returns dimension of lmim object (dimension of field A)
%
% Note the difference to size(p,dim), which is the default array size! 

% if isempty(p)
%     di = 0; 
%     nc = 0;
% else
p=lmim(p);
if nargin<2
    dim=[1 2];
end
if nargout==0
    size(p.A,dim)
end
if nargout==1
    di=size(p.A,dim);
end
if nargout==2
    [di,nc]=size(p.A,dim);
end
end
% end