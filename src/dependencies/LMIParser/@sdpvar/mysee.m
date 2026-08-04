function mysee(lmi);
%display yalmip variable

v=getvariables(lmi);
s0=getbasematrix(lmi,0);
sv=s0*0;
l=length(v);
for j=1:l
    sv=sv+getbasematrix(lmi,v(j))*v(j);
end;
if l>0
    Involved_Variables=v
    Constant_part=full(s0)
    Variable_part=full(sv)
else
    disp(' ')
    disp(full(s0))
end




