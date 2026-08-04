function s=horzcat(varargin)
% function s=horzcat(varargin)
%
% Horizontal concatenation [s1 s2 ...] of lmibl objects s1,s2 ...


%find indices of non-empty lmibl objects
iv = find(~cellfun(@isempty, varargin));

%select only nonempty entries
sinp=varargin(iv);

%check whether all are lmibl objects.
ind = sum(~cellfun(@(x) isa(x,'lmibl'), sinp));
if ind>0
    error('Can only horizontally concatenate a set of lmibl objects.')
end

n=length(sinp);
s=sinp{1};
for j=2:n
    %elements in previous so are consistent
    so=s;
    %elments in new part sn need to be checked for consistency with so
    sn=sinp{j};
    %since both so and sn are arrays of lmibl, they are consistent by
    %themselves.

    %find intersection sno of sn.na with so.na
    [sno,isn,iso]=intersect({sn.na},{so.na});

    %if sno not empty:
    %  isn are indices of sno in sn
    %  iso are FIRST indices of sno so.na.

    for k=1:numel(isn)
        %check s1 and s2
        s1=so(iso(k));
        s2=sn(isn(k));

        s1r=strcmp(s1.ty,'rep');  %yes if s1 is repeated
        s2r=strcmp(s2.ty,'rep');  %yes if s2 is repeated
        s1s=isequal(s1.di,[1 1]); %yes if s1 is sclar
        s2s=isequal(s2.di,[1 1]); %yes if s2 is sclar

        %Logical table rep (reptead) sca (scalar) for s1 versus s2
        %Concatenatio allowed (+) and not allowed (-)
        %
        %                 s2: ga         ga I        ga          [ga ga]
        %s1                   (rep,sca)  (rep,~sca)  (~rep,sca)  (~rep,~sca)
        %ga      (rep,sca)    +          +           +           -
        %ga I    (rep,~sca)   +          +           +           -
        %ga      (~rep,sca)   +          +           +           -
        %[ga;ga] (~rep,~sca)  -          -           -           check next

        %Error: s1 is not repeated, not scalar XOR s2 is not repeated, not scalar
        if xor( ~s1r && ~s1s , ~s2r && ~s2s )
            error('Variables with identical names: Repeated or scalar not matching with other.')
            %         %if s2 is not repeated and not scalar then s1 cannot be concatenated
            %         if ~s2r && ~s2s
            %             error('Variables with identical names: Repeated or scalar not matching with non-repeated and non-scalar.')
            %         end
        else
            %check next
            if and( ~s1r && ~s1s , ~s2r && ~s2s )
                %if
                if  ~strcmp(s1.ty,s2.ty)
                    error('Variables with identical names: Both non-scalar and non-repeated, TYPES do not match.')
                else
                    %check that transpose flags are equal and dimensions are equal
                    lo1= isequal(s1.tr,s2.tr) & isequal(s1.di,s2.di);

                    %check that transpose flags are opposite and swapped dimensions are equal
                    lo2= isequal(s1.tr,~s2.tr) & isequal(s1.di([2,1]),s2.di);

                    %both condition are not true, concatenation is prohibited.
                    if (~lo1) & (~lo2)
                        error('Variables with identical names: Both non-scalar and non-repeated, DIMESIONS do not match.')
                    end
                end
            end
        end
    end
    %call original matlab function
    s=builtin('horzcat', so,sn);
end

