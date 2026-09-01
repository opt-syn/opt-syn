classdef lmibl
    % Objects lmibl (LMI blocks) define basic matrix variables.
    %
    % It is discouraged to work with these objects directly. 
    % Rather generate corresponding identity maps with lmim.         
    % ---------------------------------------------------------------------  
    % function s=lmibl(na,nr,nc,type)
    %
    % Generates variable with
    % name na, dimension [nr,nc], type type.

    properties (SetAccess = private)
        % name of matrix variable
        na =''
        % type of matrix variable: {'full','sym','rep'}
        ty =''
        % Dimension [#rows,#columns] of matrix variable
        di =[]
        % Flag for transposed (logical yes=1) or not (logical 0=no)
        tr =[]
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = private)
        % Coupled validation  ---------------------------------------------
        function validate(p)
            %symmetric or repeated blocks must be square 
            if ismember(p.ty,{'sym','rep'})
                if p.di(1)~=p.di(2)
                    error('Symmetric or repeated lmibl objects must be square.');
                end
            end
        end
        %------------------------------------------------------------------
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        % Checking types if setting values --------------------------------
        function p = set.na (p, val)
            if ~isempty(val)
                if ~ischar(val)
                    error('Name must be character array.')
                end
            end
            p.na = val;
        end

        function p = set.ty (p, val)
            if ~isempty(val)
                if ~ischar(val)
                    error('Type must be character array.')
                end
                if ~ismember(val,{'full','sym','rep'})
                    error('Type not implemented.')
                end
            end
            p.ty = val;
            %validate whether dimension is consistent with type.
            p.validate;
        end

        function p = set.di (p, val)
            if ~isempty(val)
                mustBeInteger(val);
                mustBePositive(val);
                [nr,nc]=size(val);
                if nr~=1 | nc~=2
                    error('Dimension must be a positive integer row vector with 2 entries.');
                end
            end
            p.di = val;
            %validate whether dimension is consistent with type.
            p.validate;
        end

        function p = set.tr (p, val)
            if ~isempty(val)
                if ~(islogical(val) && isequal(size(val), [1 1]))
                  error('Must be a scalar logical.');
                end                
            end
            p.tr = val;
        end

        %----------------------------------------------------------------------
        function s=lmibl(na,nr,nc,type)            
            if nargin>=1                
                %default dimension 1x1 and 'full' if only name is given                
                s.na=na;
                s.di=[1 1];
                s.ty='full';
                s.tr=false;
                if nargin>=2;
                    %default square if only first dimension is given
                    s.di=[nr nr];
                end
                if nargin>=3
                    %repeated if nc is zero
                    if nc==0 & nargin<4
                        nc=nr;
                        s.ty='rep';
                    end
                    %otherwise define dimension as in input
                    s.di=[nr nc];
                end
                %take type as provided
                if nargin>=4
                    s.ty=type;
                end                
                s.validate;
            else
                %prevents from generating s.na='',s.ty='',s.di=[], s.tr=[].
                %empty lmibl are not used in the current package.
                s=lmibl.empty(0,0);
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end