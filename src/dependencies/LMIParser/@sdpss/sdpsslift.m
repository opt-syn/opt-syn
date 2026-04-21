function Gl=sdpsslift(G,N);
%Generates lifted system representation for sdpss systems
%G: original sdpss system
%N: length of lifting 
%Example N=2: [A^2] [AB B] [C;CA] [D 0;CB D]

G=sdpss(G);
Gl=sdpss;

[A,B,C,D]=sdpssdata(G);
if isempty(A)
    if isnumeric(D)
        [ny, nu] =size(D);
    else
        [ny, nu] =dim(D);
    end
    Al = [];
    Bl = zeros(0, nu*N);
    Cl = zeros(ny*N, 0);
    Dl = drep(D,N);

else
    Al=A;Bl=B;Cl=C;Dl=D;
    for j=2:N;
    
        
        Dl=[Dl zeros(dim(Dl,1),dim(D,2))];
        Dl=[Dl;C*Bl D];
        Al=A*Al;
        Bl=[A*Bl B];
        Cl=[C;Cl*A];    
    end;
end
Gl.A=Al;
Gl.B=Bl;
Gl.C=Cl;
Gl.D=Dl;
end