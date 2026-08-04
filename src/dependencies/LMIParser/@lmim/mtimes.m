function s=mtimes(s1,s2)
% function s=mtimes(s1,s2)
%
% Product of lmim objects s1 and s2
% 
% - s1 or s2 scalar map, scalar variables
% - s2 or s1 constant map of any dimension 
%
% - s1,s2 matching dimension of product 
%   s1.C*s2.B must be zero to enforce that s1*s2 is affine.   
%   Column partion of s1 and row partion of s2 may not match

%(A1 + B1*X1*C1)*(A2 + B2*X2*C2)=
%     A1*A2 + A1*B2 *X2* C2 + B1 *X1* C1*A2 + B1 *X1* C1*B2 *X2* C2

s1=lmim(s1);
s2=lmim(s2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Take care of scalar constant factors 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if s1 is constant scalar, generate s1*eye(rowdim s2)
if isempty(s1.bl) & isequal(dim(s1),[1,1])
    d=dim(s2,1);
    s1.A=s1.A*eye(d);
    s1.B=zeros(d,0);
    s1.C=zeros(0,d);
    s1.rpar=d;
    s1.cpar=d;
end
%if s2 is constant scalar, generate s2*eye(coldim s2)
if isempty(s2.bl) & isequal(dim(s2),[1,1])
    d=dim(s1,2);
    s2.A=s2.A*eye(d);
    s2.B=zeros(d,0);
    s2.C=zeros(0,d);
    s2.rpar=d;
    s2.cpar=d;
end

s=lmim;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Improve error checking: 
%   If map is scalar but variables are not scalar, 
%   current error reports mismatch of dimension. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Case s1*s2: s1=s1.A constant s2 = s2.a + s2.b*s2.bl*s2.c (all s2.bl dimension [1 1])
% A(a+bxc)=Aa+(Ab)(Ix)(Ic), I has column dimension of A

% column dimension of s1.A
d=size(s1.A,2);
%matrix with stacked column dimensions of s2.bl
if ~isempty(s2.bl)
    div=vertcat(s2.bl(:).di);
else
    div=[];
end
if d>1 && isempty(s1.bl) && isequal(size(s2.A),[1 1]) && isequal(div,ones(size(div,1),2))    
    s.A=kron(s1.A,s2.A);
    s.B=kron(s1.A,s2.B);
    s.C=kron(eye(d),s2.C);
    s.rpar=size(s1.A,1);
    s.cpar=d;
    s.bl=lmibl(s2.bl(1).na,d,0);
    for j=2:length(s2.bl);
        s.bl(j)=lmibl(s2.bl(j).na,d,0);
    end
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Case s1*s2: s1 = s1.a + s1.b*s1.bl*s1.c (all s2.bl dimension [1 1]) and s2=s2.A constant  
    % (a+bxc)A=aA+(bI)(xI)(cA), I has row dimension of A    
    
    % row dimension of s2.A
    d=size(s2.A,1);
    %matrix with stacked column dimensions of s1.bl
    if ~isempty(s1.bl)
        div=vertcat(s1.bl(:).di);
    else
        div=[];
    end
    if d>1 && isempty(s2.bl) && isequal(size(s1.A),[1 1]) && isequal(div,ones(size(div,1),2))
        s.A=kron(s1.A,s2.A);
        s.B=kron(s1.B,eye(d));
        s.C=kron(s1.C,s2.A);
        s.rpar=d;
        s.cpar=size(s2.A,2);        
        s.bl=lmibl(s1.bl(1).na,d,0);
        for j=2:length(s1.bl);
            s.bl(j)=lmibl(s1.bl(j).na,d,0);
        end
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %All other cases
        [k1,m1]=dim(s1);
        [k2,m2]=dim(s2);
        if m1~=k2;
            error('Dimensions not compatible.')
        else
            if norm(s1.C*s2.B)~=0
                error('Product not affine.')
            else
                % concatenation of two block cells is diagonal combination
                s=s1;
                s.A=s1.A*s2.A;
                s.B=[s1.B s1.A*s2.B];
                s.C=[s1.C*s2.A;s2.C];
            end
            %row partition of s1 and column partition of s2 inherited
            s.rpar=s1.rpar;
            s.cpar=s2.cpar;
            %put at end for validation
            s.bl=[s1.bl s2.bl];
        end
    end
end

