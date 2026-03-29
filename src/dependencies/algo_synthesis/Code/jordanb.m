function j = jordan(n)
%determines upper Jordan block of dimension n with eigenvalue zero
if n==0;
    j=[];
elseif n==1;
    j=0;
else
   n=n-1;
   j=[zeros(n,1) eye(n);0 zeros(1,n)];
end;