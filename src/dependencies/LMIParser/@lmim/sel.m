function s=sel(p,indo,indi)
%s=sel(p,indo,indi)
%
%Select channels from lmim object indexing.
%Partitioning is lost for simplicity if selection takes place. 
%Command s=sel(p) removes any partitioning information. 
%
%Note that variable structure remains untouched!
%Variable adaptation (dimensionality reduction) not implemented.
%
%Example s=sel(p,[1 3 2],'1:end-3');

p=lmim(p);
s=p;
if nargin>1
    s.A=p.A(indo,indi);
    s.B=p.B(indo,:);
    s.C=p.C(:,indi);
    %warning('Partitions not interited.')
end
s.rpar=size(s.A,1);
s.cpar=size(s.A,2);
end