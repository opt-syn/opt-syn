function sout=sel(s,indo,indi)
%sout=sel(s,indo,indi)
%
%Select channels from sdpss system s by indexing.
%Numeric indexing for indo and inti standard.
%Other indexing in matlab arrays used as strings.
%
%Example sout=sel(s,[1 3 2],'1:end-3');

s=sdpss(s);
if nargin<=1
    sout=s;
else
    [k,m]=size(s.D);
    Io=eye(k);
    Ii=eye(m);

    so=sdpss;
    si=sdpss;
    si.D=Ii;
    if isstr(indo);
        eval(['so.D=Io(' indo ',:);']);
    else
        so.D=Io(indo,:);
    end;
    if nargin>2;
        if isstr(indi);
            eval(['si.D=Ii(:,' indi ');']);
        else
            si.D=Ii(:,indi);
        end;
        sout=so*s*si;
    end
end