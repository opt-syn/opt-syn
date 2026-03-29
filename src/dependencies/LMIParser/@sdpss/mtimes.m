function s=mtimes(s1,s2)
%Series interconnection (product) of two sdpss objects s1 and s2
s1=sdpss(s1);
s2=sdpss(s2);

s=sdpss;
n1=size(s1.A,1);
n2=size(s2.A,1);
if n1==0;
    s.A=s2.A;
    s.B=s2.B;
    
    if isempty(s2.C)
        s.C = s2.C;
    else
        s.C=s1.D*s2.C;
    end
    s.D=s1.D*s2.D;
else
    if n2==0;
        s.A=s1.A;
        s.B=s1.B*s2.D;
        s.C=s1.C;
        s.D=s1.D*s2.D;
    else
        s.A=[s1.A s1.B*s2.C;zeros(n2,n1) s2.A];
        s.B=[s1.B*s2.D;s2.B];
        s.C=[s1.C s1.D*s2.C];
        s.D=s1.D*s2.D;
    end
end
end
