function s=parcomb(par1,par2)
% function s=parcomb(par1,par2)
%
% Cobmine partition vectors (containing dimensions) by coarsening:
% par1=[d1 d2 ... dN]
% par2=[e1 e2 ... eM]

arguments
    par1 double {mustBeInteger, mustBeNonnegative, mustBeVector(par1,"allow-all-empties")}
    par2 double {mustBeInteger, mustBeNonnegative, mustBeVector(par2,"allow-all-empties")}
end

if sum(par1)==0    
    par1=par2;
else
    if sum(par2)==0;
        par2=par1;
    else
        if sum(par1)~=sum(par2)
            error('Dimensions corresponding to partitions not equal.')
        end

        ind=1;
        while ind<=length(par1) & ind<=length(par2)
            while par1(ind)~=par2(ind)
                %if d(ind)<e(ind) then combine d(ind) and d(ind+1) by addition
                %shortens vector d
                if par1(ind)<par2(ind) & ind<length(par1)
                    par1(ind)=par1(ind)+par1(ind+1);
                    par1=par1([1:ind,ind+2:end]);
                end
                %if d(ind)>e(ind) then combine e(ind) and e(ind+1) by addition
                %shortens vector e
                if par1(ind)>par2(ind) & ind<length(par2)
                    par2(ind)=par2(ind)+par2(ind+1);
                    par2=par2([1:ind,ind+2:end]);
                end
                %repeat this step until d(ind)=e(ind)
                %this ends up with lists of one entry containing the dimension
            end
            ind=ind+1;
        end
    end
end

if isequal(par1,par2)
    s=par1;
else
    %this should never happen after dimension check at beginning
    error('Combination of partitions not possible.')
end