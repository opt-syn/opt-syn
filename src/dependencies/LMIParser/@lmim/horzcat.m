function s=horzcat(varargin)
% function s=horzcat(varargin)
%
% Horizontal concatenation [s1 s2 ...] of VALUES of lmim objects s1,s2 ...
% Row partition is coarsening of those in s1,s2,...
%
% Note the difference to concatenation of lmim arrays! 


n=length(varargin);
sinp=varargin;
if isempty(sinp{1})
    s = [];
else
    s=lmim(sinp{1});
end
for j=2:n
    s1=s;
    if ~isempty(sinp{j})
        s2=lmim(sinp{j});
        if isempty(s)
            s = s2;
        else            
            if size(s1.A,1)~=size(s2.A,1)
                error('Row dimensions of maps not all equal.')
            end
            s.A=[s1.A s2.A];
            s.B=[s1.B s2.B];
            s.C=blkdiag(s1.C,s2.C);
            s.cpar=[s1.cpar s2.cpar];
            %combine row partion
            s.rpar=parcomb(s1.rpar,s2.rpar);
            %put at end for validation
            s.bl=[s1.bl s2.bl];
        end
    end
end
end

