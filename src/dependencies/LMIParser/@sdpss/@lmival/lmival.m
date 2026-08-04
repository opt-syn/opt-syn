classdef lmival
    % Objects lmival are row vectors of doubles.
    % Accessing follows the default behavior of matlab.
    % Using the constructor reshapes arrays to rows. 
    % Only directly assing with ONE index, like obj(4)=rand(3,2)
    % -----------------------------------------------------------------------
    % function s = lmival(p,val)
    %
    % Only p: Generates lmival s with value p (double, cell of doubles, lmival array)
    % Otherwise append val (double, cell of doubles, lmival array) to lmival array p
        
    properties
        %value of each entry in vector 
        val double = [];
    end

    % Pre-computed constants or internal states
    properties (Access = private)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Static,Access = private)
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        %------------------------------------------------------------------
        function p=set.val(p,val)            
            p.val=val;                        
        end
        %------------------------------------------------------------------
        function s=ctranspose(p)
            %not allowed
            error('No transposition of lmival input.')
        end
        function s=transpose(p)
            %not allowed
            error('No transposition of lmival input.')
        end
        function s=vertcat(p1,p2)
            %not allowed
            error('No vertical concatenation of lmival inputs.')
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function s = lmival(p,val)
            % function s = lmival(p,val)
            %
            % Only p: Generates lmival s with value p (double, cell of doubles, lmival array)
            % Otherwise append val (double, cell of doubles, lmival array) to lmival array p
            
            %translate input into lmival array
            if nargin==1
                s=checkval(p);
            end
            %append second input to first input which must be lmival array
            if nargin==2;
                if ~isa(p,'lmival')
                    error('First argument must be lmival.')
                end
                %empty val not allowed (could be changed to no action)
                if isempty(val)
                    error('Second argument should not be empty.')
                end
                %translate val into lmival array
                sh=checkval(val);
                le=numel(sh);
                s=p;
                %append to row vector
                s=reshape(p,1,[]);
                s(end+(1:le))=sh;                
            end            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        function s=double(p)
            % function s=double(p)  
            % Generates single double or cell array of doubles
            if numel(p)==1
                s=p.val;
            else                
                for j=1:numel(p)
                    s{j}=p(j).val;
                end
            end
        end
        %------------------------------------------------------------------
        function s=eig(p)
            % function s=eig(p)
            % Generates lmival for eigenvalues in lmival list p (if square matrices)
            s=lmival;
            for j=1:numel(p)
                A=double(p(j));
                if size(A,1)==size(A,2)
                    s(j)=eig(A);
                else
                    s(j)=lmival([]);
                end
            end
        end
        %------------------------------------------------------------------
        function disp(p)
            disp(['List of ' num2str(numel(p)) ' of LMI value(s).'])
            disp(' ')
            for j=1:numel(p);
                if isempty(p(j).val)
                    disp('    []');
                    disp(' ');
                else
                    disp(p(j).val)
                end
            end
        end
        %------------------------------------------------------------------
        function s=dim(p)                        
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end