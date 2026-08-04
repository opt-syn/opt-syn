function so=checkval(val)
% function so=checkval(val)
%
% val (double, cell of doubles, @lmival array) is checked for integrity and
% changed into lmival object array
if isempty(val)
    so=lmival;
else
    %empty object
    so=lmival;
    switch class(val)
        %double just assigned to val
        case 'double'
            so.val=val;
        %reshape cell array into cell row
        %collect into lmival array so
        case 'cell'
            val=reshape(val,1,[]);
            so(1)=val{1};
            for j=2:numel(val)
                so(end+1)=val{j};
            end
        %reshape lmival array into lmival row 
        case 'lmival'
            val=reshape(val,1,[]);
            so=val;
        otherwise
            error('Type of input values not supported.')
    end
end

