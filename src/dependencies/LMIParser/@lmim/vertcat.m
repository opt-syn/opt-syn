function s=vertcat(varargin)
% function s=vertcat(varargin)
%
% Vertical concatenation [s1;s2;...] of VALUES of lmim objects s1,s2,....
% Column partition lost is coarsening of those in s1,s2,...
%
% Note the difference to concatenation of lmim arrays! 

n=length(varargin);
sinp=varargin;
s=lmim(sinp{1});
for j=2:n
    s1=s;
    s2=lmim(sinp{j});
    if size(s1.A,2)~=size(s2.A,2)
        error('Column dimensions of maps not all equal.')
    end
    s.A=[s1.A;s2.A];
    s.B=blkdiag(s1.B,s2.B);
    s.C=[s1.C;s2.C];    
    s.rpar=[s1.rpar s2.rpar];
    %combine column partion
    s.cpar=parcomb(s1.cpar,s2.cpar);
    %put at end for validation 
    s.bl=[s1.bl s2.bl];
end
end