function [so,info]=yalsolve(lmis,opt)
% function [so,info]=yalsolve(lmis,opt)
%
% Defines lmi for yalmip to solve SDP.
% lmis is system of LMIs object
% opt is option according to sdpsettings.

if ~isa(lmis,'lmis')
    error('First argument must be an LMI sytem object lmis.')
end

if nargin<2
    opt=sdpsettings;    
    opt.solver='sdpt3';
    opt.verbose=0;
end

clear('yalmip')
lmisys=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define variables
for iv=1:length(lmis.var);
    %iv index of variable in varl(iv) and vas{iv}
    na=char(lmis.var(iv));
    ty=lmis.bl(iv).ty;
    if strcmp(ty,'rep')
        di=num2str(lmis.bl(iv).di(1));
        eval([na '= sdpvar(1,1,''full'');']);
        eval([na '= ' na '*eye (' num2str(di) ');']);
    else
        nr=num2str(lmis.bl(iv).di(1));
        nc=num2str(lmis.bl(iv).di(2));
        eval([na '= sdpvar(' nr ',' nc ',' ''''  ty '''' ');']);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set up lmi system for yalmip: il is index of LMI.
%Ignore partitions. Should be fine.

for il=1:length(lmis.lmim);
    s=lmis.lmim(il);
    lm=s.A;
    Bi=0;
    Ci=0;
    for j=1:length(s.bl);
        na=s.bl(j).na;
        %check whether in varl: this should be always true!
        [ind,iv]=ismember(na,lmis.var);
        if ind==0
            error('Variable list lmis.var not properly constructed.')
        else
            nr=s.bl(j).di(1);
            nc=s.bl(j).di(2);
            L=s.B(:,Bi+(1:nr));
            R=s.C(Ci+(1:nc),:);
            Bi=Bi+nr;
            Ci=Ci+nc;
            if s.bl(j).tr
                %actual variable is tranposed version of implemented one
                %implement trasposed variable
                na=[ na ''''];
            end
            eval(['lm = lm + L*' na '*R;']);
        end
    end
    %use outer factor if exisiting
    R=double(lmis.fac(il));
    if isempty(R)
        nc=dim(s,2);
        R=eye(nc);        
    end    
    lmisys=[lmisys,R'*lm*R<=0];
end
nvar=yalmip('nvars');
%disp(['Number of decision variables: ' num2str(nvar)])
if nvar>10000;
    warning('on')
    warning('More that 10000 variables.')
end

if isempty(lmis.cost)
    co=[];
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %set up cost vector
    s=lmis.cost;
    %cost is always 1x1 by default.
    co=s.A;
    Bi=0;
    Ci=0;
    for j=1:length(s.bl)
        %name of current variable        
        na=s.bl(j).na;
        [ind,iv]=ismember(na,lmis.var);
        if ind==0
            error('Variable list lmis.var not properly constructed.')
        else
            nr=s.bl(j).di(1);
            nc=s.bl(j).di(2);
            L=s.B(:,Bi+(1:nr));
            R=s.C(Ci+(1:nc),:);
            Bi=Bi+nr;
            Ci=Ci+nc;
            if s.bl(j).tr
                %actual variable is tranposed version of implemented one
                %implement transposed variable                
                na=[ na ''''];
            end
            eval(['co = co + L*' na '*R;']);
        end
    end
end

tic
dia=optimize(lmisys,co,opt);
info.time=toc;

so=lmis;
%return info about specific implementation in yalmip
info.lmisys=lmisys;
info.decnbr=nvar;
info.lhs=check(lmisys);
info.lhsmax=-max(info.lhs);
info.dia=dia;
if info.lhsmax>0
    warning('on')
    warning('LMIs seem not feasible. Double check solution!')
end
%return solution in lmis structure
so.val=lmival(so.val);
for iv=1:length(lmis.var)
    eval(['so.val(iv)=double( ' char(lmis.var(iv)) ');']);
end
so.dia=[double(co) info.lhsmax];
end

