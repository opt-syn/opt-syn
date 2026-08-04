function s=plus(s1,s2)
%Parallel interconnection (sum) of two sdpss objects s1 and s2
s1=sdpss(s1);
s2=sdpss(s2);
s=sdpss;
n1=size(s1.A,1);
n2=size(s2.A,1);
s.A=[s1.A zeros(n1,n2);zeros(n2,n1) s2.A];
s.B=[s1.B;s2.B];
s.C=[s1.C s2.C];
s.D=s1.D+s2.D;
end
