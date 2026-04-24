function [so,info]=lmisolve(lmis,opt)
% function [so,info]=lmisolve(lmis,opt)
%
% Defines lmi for LMIlab (Matlab) and solves it
% lmis is system of LMIs object
% opt is option according LMIlab.

%lmis is system of LMIs object
%defines lmi for matlab solver and solves it

% warning('on')
warning('off', 'Robust:lmi:NonSymmetricTerm1')


if ~isa(lmis,'lmis')
    error('First argument must be an LMI sytem object lmis.')
else
    check(lmis)
end


if nargin<2
    %opt=[0 0 0 0 1];
    reltol=1e-4;
    maxiter=300;
    feasradius=-1; %flexible bound method
    L=20;
    verbose=1;
    opt=[reltol maxiter feasradius L verbose];
end

setlmis([]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define variables
%vai(iv) should be iv
%vas{iv} collects structure in definition of LMI system
vai=[];
vas={};
vas=lmival;
for iv=1:numel(lmis.var);
    %iv index of variable in varl(iv) and vas(iv)
    nr=lmis.bl(iv).di(1);
    nc=lmis.bl(iv).di(2);
    switch lmis.bl(iv).ty
        case 'full'
            [vai(iv),~,vas(iv)]=lmivar(2, [nr nc]);
        case 'sym'
            [vai(iv),~,vas(iv)]=lmivar(1,[nr 1]); %full symmetric, guaranteed to be square by construction
        case 'rep'
            [vai(iv),~,vas(iv)]=lmivar(1,[nr 0]); %nr times repeated, guaranteed to be square by construction
        otherwise
            error('Type not implemented.')
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set up lmi system: il is index of LMI
nterms=0;
for il=1:numel(lmis.lmim);
    s=lmis.lmim(il);
    %extract list of coefficients in actual LMI map s
    [A,B,C]=coelist(s);
    for i=1:numel(s.cpar)
        for k=i:numel(s.cpar)
            if norm(A{i,k})>0
                %implement constant term
                lmiterm([il,i,k,0],A{i,k});
                nterms=nterms+1;
            end
            %list of indices for variable selection
            %this list is turned empty by implmenting blocks
            %can implement symmetrized blocks in this way
            ivarl=1:numel(s.bl);

            %only proceed if list not empty
            while numel(ivarl)>0
                j=ivarl(1);

                %first entry is actual variable index j
                %remove first entry from list
                ivarl=ivarl(2:end);

                %actual left and right factor in map
                L=B{i,j};
                R=C{j,k};

                %if both zero, nothing implemented and only ivarl reduced
                % if norm(L)*norm(R)~=0
                if any(L, 'all') && any(R, 'all')
                    %extract name of current variable
                    na=s.bl(j).na;

                    %check whether in varl: this should be always true!
                    [ind,iv]=ismember(na,lmis.var);
                    if ind==0
                        error('Variable list lmis.var not properly constructed.')
                    end

                    %if current block transposed define factor=-1 to adjust
                    %the sign of vai(iva) in block implementation

                    factor=1;
                    if s.bl(j).tr
                        %current block transposed
                        factor=-1;
                    end

                    %flag for symmetrization
                    lo=false;

                    %search for symmetric counterpart (on diagonal only)
                    if i==k
                        %run through remainin list
                        for rest=1:numel(ivarl)
                            jnew=ivarl(rest);
                            Lnew=B{i,jnew};
                            Rnew=C{jnew,k};

                            %has it the same name as actual variable?
                            lo1=strcmp(na,s.bl(jnew).na);

                            %is actual left factor transpose of new right factor?
                            lo2=isequal(L,Rnew');

                            %is actual right factor transpose of new left factor?
                            lo3=isequal(R,Lnew');

                            %is actual variable transpose of new variable?
                            %is s.bl(j) transpose of s.bl(jnew)?
                            lo4=false;
                            if (s.bl(j).tr==not(s.bl(jnew).tr));
                                lo4=true;
                            else
                                %even if the flags .tr are equal, lo4 still
                                %true if both variables are of type symmetric
                                %or repeated
                                if strcmp(s.bl(jnew).ty,'sym') | strcmp(s.bl(jnew).ty,'rep');
                                    lo4=true;
                                end
                            end

                            %symmetric counterpart found: lo=true
                            lo=lo1 & lo2 & lo3 & lo4;

                            %if lo true remove index jnew from list to
                            %avoid later implementation of symmetric
                            %counterpart
                            if lo
                                ivarl=setdiff(ivarl,jnew);
                                %break loop if symmetric counterpart found
                                break
                            end
                        end
                    end

                    %if L is identity then set L to 1
                    if isequal(L,eye(size(L,1)))
                        L=1;
                    end
                    %if R is identity then set R to 1
                    if isequal(R,eye(size(R,2)))
                        R=1;
                    end

                    %lo can only be true for diagonal blocks (i=k)!
                    if lo
                        %if lo is true then use symmetric implementation
                        lmiterm([il,i,k,factor*vai(iv)],L,R,'s');
                    else
                        %if lo is false then no symmetric implementation
                        lmiterm([il,i,k,factor*vai(iv)],L,R);
                    end
                    nterms=nterms+1;
                end
            end
        end
    end
    %set outer factors if exisiting
    R=double(lmis.fac(il));
    if ~isempty(R)
        lmiterm([il 0 0 0],R);
    end
end
lmisys=getlmis;
nvar=decnbr(lmisys);
%disp(['Number of decision variables: ' num2str(nvar)])
%disp(['Number of terms: ' num2str(nterms)])
if nvar>1000;
    warning('on')
    warning('More that 1000 variables.')
end

if isempty(lmis.cost)
    tic
    [vopt,xopt]=feasp(lmisys,opt);
    info.time=toc;
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %set up cost vector by handling constant part of cost map directly
    %recall that cost map is guaranteed to be 1x1
    %c is vector for LMIlab implementation
    c=zeros(nvar,1)';
    s=lmis.cost;
    [A,B,C]=coelist(s);
    for k=1:nvar
        %run over variables in cost map 
        for j=1:numel(s.bl)
            %name of current variable
            na=s.bl(j).na;
            [ind,iv]=ismember(na,lmis.var);
            if ind==0
                %this should never occur if lmis constructed properly.
                error('Variable list lmis.var not properly constructed.')
            else
                %if variable in list then the index of the variable in the 
                %LMI implmentation is also iv by constrution. We can
                %extract the corresponding k-basis matrix Xv with defcx.
                %
                %requires transposition if actual variable is transposed
                Xv=defcx(lmisys,k,iv);
                if s.bl(j).tr
                    Xv=Xv';
                end
                c(k)=c(k)+B{1,j}*Xv*C{j,1};
            end
        end
    end
    tic
    [vopt,xopt]=mincx(lmisys,c,opt);
    info.time=toc;
end

so=lmis;
so.status = isempty(vopt);
if ~so.status;
    %return info about specfic implementation in LMILab
    info.opt=opt;
    info.lmisys=lmisys;
    info.xopt=xopt;    
    info.vopt=vopt;
    info.decnbr=nvar;
    info.nterms=nterms;
    info.lhsmax=[];
    if ~isempty(lmis.cost)
        info.c=c;
    end
    
    lmivarl=evallmi(lmisys,xopt);
    for il=numel(lmis.lmim)
        [lhs,rhs]=showlmi(lmivarl,il);
        info.lhsmax=[info.lhsmax max(eig(lhs))];
    end
    info.vas=vas;
    if max(info.lhsmax)>0
        so.status = 1;
        % warning('on')
        % warning('LMIs seem not feasible. Double check solution!')
    end
    
    
    %return solution in lmis structure
    so.val=lmival(so.val);
    for iv=1:numel(lmis.var)
        so.val(iv)=dec2mat(lmisys,xopt,vai(iv));
    end
    %add constant term of cost, ignored so far
    if isempty(lmis.cost)
        so.dia=vopt;
    else
        so.dia=[lmis.cost.A+vopt info.lhsmax];
    end
end
