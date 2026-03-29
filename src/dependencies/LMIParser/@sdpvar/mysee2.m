function mysee(lmi);
%display yalmip variable 

v=getvariables(lmi);
s0=getbasematrix(lmi,0);
sv=s0*0;
for j=1:length(v);
    sv=sv+getbasematrix(lmi,v(j))*v(j);
end;
Involved_Variables=v
Constant_part=full(s0)
Variable_part=full(sv)



