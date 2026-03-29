    classdef lmim
    % Objects lmim are affine maps to build lmis
    %
    % A+B*X*C   A,B,C describing matrices and X is block diagonal variable
    %           with block information in array bl of lmibl objects 
    % A         Constant map: X empty, B has zero columns, C has zero rows
    %----------------------------------------------------------------------
    % s=lmim(p,rpar,cpar)
    %
    % p is lmim:   Only change partitions to rpar, cpar
    %              No array of lmim allowed. 
    %              This is the only way to adapt the partiton.
    % p is lmibl:  Generate indentity map (with partition)
    %              No array of lmibl allowed. 
    % p is double: Generate constant map (with partition)           
    %              No array of double allowed. 
    % p is char:   Construct lmibl object from input and 
    %              proceed as if p is lmibl.

    properties
        %constant coefficient matrix
        A = []
        %left coefficient 
        B = []
        %right coefficient 
        C = []
        %array of lmibl objects (LMI variables on diagonal)
        bl  = []
        %row partition [dim1 dim2 ... dimN] of map for implementation in LMIlab               
        rpar = []
        %column partition [dim1 dim2 ... dimN] of map for implementation in LMIlab        
        cpar = []
    end

    % Pre-computed constants or internal states
    properties (Access = private)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        % Checking types if setting values --------------------------------
        function p = set.bl (p, val)
            if ~isempty(val)
                if ~isa(val,'lmibl')
                    error('Field bl must be an lmibl object.')
                end
            end
            p.bl = val;
            p.check;
        end
        %------------------------------------------------------------------
        function obj = set.rpar(obj, val)
            if ~isempty(val)
                mustBeInteger(val);
                mustBePositive(val);
                [nr,nc]=size(val);
                if nr~=1
                    error('Field rpar must be a row vector');
                end
            end
            obj.rpar = val;
        end
        %------------------------------------------------------------------
        function obj = set.cpar(obj, val)
            if ~isempty(val)
                mustBeInteger(val);
                mustBePositive(val);
                [nr,nc]=size(val);
                if nr~=1
                    error('Field cpar must be a row vector');
                end
            end
            obj.cpar = val;
        end
        %------------------------------------------------------------------
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = protected)
        % Coupled validation ----------------------------------------------
        function p = check(p)            
            [k,m]=size(p.A);
            %check row partition 
            if k~=sum(p.rpar)
                error('Row partition of map not compatible with dimension.')
            end
            %check column partition 
            if m~=sum(p.cpar)
                error('Column partition of map compatible with dimension.')
            end      
            %check compatibility of dimensions in product B*X*C
            [nr,nc]=sizebl(p);            
            if nr~=size(p.B,2) | nc~=size(p.C,1)
                error('Inner dimension in A+B*X*C (products) not compatible.');
            end
        end
        function [nr,nc]=sizebl(p)
            nr=0;nc=0;
            for j=1:length(p.bl)
                nr=nr+p.bl(j).di(1);
                nc=nc+p.bl(j).di(2);
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
    methods
        %------------------------------------------------------------------
        function s=lmim(varargin)
            %s=lmim(lmim,rpar,cpar)
            %   s=lmim
            %   BOTH rpar,cpar exisiting: s.rpar=rpar,s.cpar=cpar
            %   This is the only way to adjust partitions!
            %s=lmim(lmibl,...)
            %   s=identity map with variable blk
            %   Ignore other arguments
            %
            %s=lmim(double,...)
            %   s=constant map with s.A=double
            %   Ignore other arguments
            %
            %s=lmim(char,...)
            %   s=identity map with varible defined by lmibl(char,...)            %  
            
            le=numel(varargin);            
            if le==0
%               s=lmim.empty(0,0);
            else                
                p=varargin{1};                   
                
                %avoid more inputs if first is lmibl or double
                if nargin>1 
                    if isa(p,'lmibl')
                        error('To many inputs if first argument is lmibl.')
                    end
                    if isa(p,'double')                        
                    error('To many inputs if first argument is double.')
                    end
                end
                if isa(p,'lmim');
                    if numel(p)>1
                        error('First lmim argument MUST NOT be an array.')
                    end
                    s=p;
                    if le==3;
                        s.rpar=varargin{2};
                        s.cpar=varargin{3};
                    end
                    s=check(s);
                end
                if isa(p,'char')
                    p=lmibl(varargin{:});
                end                
                %first entry lmibl                
                if isa(p,'lmibl')
                    if numel(p)>1
                        error('First lmibl argument MUST NOT be an array.')
                    end
                    nr=p.di(1);
                    nc=p.di(2);
                    s.A=zeros(nr,nc);
                    s.B=eye(nr);
                    s.C=eye(nc);
                    s.rpar=nr;
                    s.cpar=nc;
                    %first entry without brackets
                    s.bl=p;
                end                
                if isa(p,'double')
                    %generate constant map with empty block of correct dimension
                    if ~isempty(p)
                        s.A=p;
                        [k,m]=size(s.A);
                        s.B=double.empty(k,0);
                        s.C=double.empty(0,m);
                        s.rpar=k;
                        s.cpar=m;
                        s.bl=[];
                    end
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function mydisp(p)
            disp([ '['  num2str(sum(p.rpar)) ','  num2str(sum(p.cpar)) ']-dimensional lmim object with partitions [' num2str(p.rpar) '], [' num2str(p.cpar) '] and variables'] )
            fprintf('\t%-10s\t%s\t%-5s\t%s\n', 'Name', 'Dimension', 'Type', 'Transposed');
            for j=1:length(p.bl);
                fprintf('\t%-10s\t[%d,%d]\t\t%-5s\t%d\n', p.bl(j).na, p.bl(j).di, p.bl(j).ty, p.bl(j).tr);
            end
        end
        function disp(pi)
            %short display of array of lmim objects
            for j=1:length(pi);
                p=pi(j);
                mydisp(p)
            end
            disp('-------------------------------------------------------------------------------')
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = protected)
        %------------------------------------------------------------------
        function so=subsloc(s,var,val)
            %evaluate s in variables given in var by values in val
            %
            %s    is lmim object (no array)
            %var  is list of variables (see function exvar) 
            %val  is list of values (see function checkval)
            %
            %so   is lmim objec in which the variables are substituted with values
            %
            %internal function, no checking of number of inputs

            %turn list of values (see checkval) into @lmival array
            val=checkval(val);
            %turn list of variables (see exvar) into @string array
            var=exvar(var);

            if numel(val)==0
                error('No values detected.')
            end
            if numel(val)~=numel(var)
                error('Number of variables must match number of values.')
            end

            %initialize
            s=lmim(s);
            so=s;

            %collect info about column and row indices of variable in s.B and s.C
            Bi=0;
            Ci=0;
            Bind=[];
            Cind=[];
            %collect info about indices of si.bl
            blind=[];

            for j=1:numel(s.bl)
                %block dimension
                nr=s.bl(j).di(1);
                nc=s.bl(j).di(2);

                %is current variable element of list var?
                [ind,iv]=ismember(s.bl(j).na,var);

                %if current variable in list then substitute
                if ind~=0
                    %pick value related to current variable
                    sval=double(val(iv));
                    if s.bl(j).tr
                        %transpose value if variable is transposed in map
                        sval=sval';
                    end
                    %update constant
                    %skip indices if variable substituted

                    %%% TYPE CHECKING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    switch s.bl(j).ty
                        case 'sym'
                            %must be square
                            if size(sval,1)~=size(sval,2)
                                error(['Value for ' '''' char(var(iv)) '''' ' not square.']')
                            else
                                %must be symmetric
                                if norm(sval'-sval)>0
                                    error(['Value for ' '''' char(var(iv)) '''' ' not symmetric.']')
                                end
                            end
                        case 'rep'
                            %must be square
                            if size(sval,1)~=size(sval,2)
                                error(['Value for ' '''' char(var(iv)) '''' ' not square.']')
                            else
                                %adapt if value is scalar by variable is not
                                if nr>1 & size(sval,1)==1 %then also nc>1 and size(sval,2)==1
                                    sval=sval*eye(nr);
                                end
                                %must be repeated times identity
                                if norm(sval-sval(1)*eye(size(sval,1)))>0
                                    error(['Value for ' '''' char(var(iv)) '''' ' not repeated.']')
                                end
                            end
                    end
                    if nr~=size(sval,1) | nc~=size(sval,2)
                        error(['Dimension of value for ' '''' char(var(iv)) '''' ' not correct.'])
                    end
                    %%% TYPE CHECKING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    so.A=so.A+s.B(:,Bi+(1:nr))*sval*s.C(Ci+(1:nc),:);

                else
                    %keep in index list if variable not substituted
                    Bind=[Bind Bi+(1:nr)];
                    Cind=[Cind Ci+(1:nc)];
                    blind=[blind j];
                end
                Bi=Bi+nr;
                Ci=Ci+nc;
            end
            so.B=s.B(:,Bind);
            so.C=s.C(Cind,:);
            so.bl=s.bl(blind);
        end        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end




