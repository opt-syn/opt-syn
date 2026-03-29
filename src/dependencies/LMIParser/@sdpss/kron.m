    function s=kron(s1,s2)
%Kronecker product of system with matrix.
s1=sdpss(s1);
s2=sdpss(s2);
s=sdpss;    
if size(s2.A,1)==0;    
    [k,m]=size(s2.D);    
    s.A=kron(s1.A,eye(k));
    s.B=kron(s1.B,s2.D);
    s.C=kron(s1.C,eye(k));
    s.D=kron(s1.D,s2.D);
else
    if size(s1.A,1)==0;
        [k,m]=size(s1.D);        
        s.A=kron(eye(m),s2.A);
        s.B=kron(eye(m),s2.B);
        s.C=kron(s1.D,s2.C);
        s.D=kron(s1.D,s2.D);
    else
        error('Both input objects have a nomempty A entry.')
    end
end
end