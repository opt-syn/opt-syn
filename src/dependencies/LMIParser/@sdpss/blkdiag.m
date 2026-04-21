function s=blkdiag(varargin)
% Block diagonal combination of sdpss objects.
n=length(varargin);
sinp=varargin;
s=sdpss;
for j=1:n
    s1=s;
    s2=sdpss(sinp{j});
    
    
    if ~isempty(s1.D)
        if isnumeric(s1.D)
            nD1 = size(s1.D,1);
        else
            nD1 = dim(s1.D,1);
        end
        s.D=blkdiag(s1.D,s2.D);
    else
        nD1 = 0;
        s.D = s2.D;
    end


    n1 = ~isempty(s1.A);
    n2 = ~isempty(s2.A);
    if n1 > 0 && n2 > 0
        s.A=blkdiag(s1.A,s2.A);
        s.B=blkdiag(s1.B,s2.B);
        s.C=blkdiag(s1.C,s2.C);
        
    elseif n1 > 0 && n2 == 0 && ~isempty(s2.D)
        %static second system
        
        if isnumeric(s2.D)
            [ny2, nu2] = size(s2.D);
        else
            [ny2, nu2] = dim(s2.D);
        end

        [~,m1]=size(s1.A);
        s.A=s1.A;
        ZB = zeros(m1, nu2);
        s.B = [s1.B, ZB];


        if isnumeric(s2.C)
            ZC = (zeros(ny2, m1));
        else
            ZC = lmim(zeros(ny2, m1));
        end
        s.C = [s1.C; ZC];



    elseif n1 == 0 && n2 > 0
        %static first system
        
        
        [~,m2]=size(s2.A);
        s.A=s2.A;
        if nD1 == 0
            s.B = s2.B;
            s.C = s2.C;
        else
            ZB = zeros(m2, nD1);
            s.B = [ZB, s2.B];
    
            if isnumeric(s2.C)
                ZC = (zeros(nD1, m2));
            else
                ZC = lmim(zeros(nD1, m2));
            end
            s.C = [ZC; s2.C];
        end
        % s.B=blkdiag(s1.B,s2.B);
        % s.C=blkdiag(s1.C,s2.C);
    end
end
end