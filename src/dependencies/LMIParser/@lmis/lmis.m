classdef lmis
    % Objects lmis are LMI systems
    %-----------------------------------------------------------------------------------------------------
    % function so=lmis(p,s,c)
    %
    % If p is lmis and c double: Adjoin (s,c) to (p.lmim,p.fac)
    %              and c char:   Set p.cost to s.  
    % If p is lmim and s double: Generate so with (so.lmim,so.fac)=(p,s)
    %                  s char:   Generate so with so.cost=s
    
    % Public, tunable properties
    properties
        % array of affine maps (lmim objects)
        lmim  = []
        % array of doubles (lmival objects): right outer factors in LMIs
        fac  = []
        % cost (lmim oject) for lmi problem
        cost  = []
        % array of strings with names of involved variables (unique names, NOT transposed)
        var = []
        % array of variables (lmibl object) corresponding to var
        bl  = []
        % value of variable in list var, as e.g. after optimization
        val = []
        %summary of diagnosis of results (value of cost and maximum of all lhs eigenvalues)
        dia = []
        %was the problem feasible (0) or infeasible (1)
        status = []
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = protected)
            % Coupled validation ----------------------------------------------
            function check(p)
                if numel(p.lmim)~=numel(p.fac)
                        error(['Number of maps and number of outer factors do not match.'])
                end
                for j=1:numel(p.lmim);
                    nc=dim(p.lmim(j),2);                    
                    [nrf,ncf]=size(p.fac(j).val);                    
                    if ncf>0 && nc~=nrf
                        error(['Entry ' num2str(j) ': Row dimension of outer factor not equal to dimension of map.'])
                    end
                end
            end
        end
    % Checking types if setting values --------------------------------
    methods
        function p = set.lmim (p, s)
            % The assignment
            %     p.lmim(end+1)=s
            % for lmis object p and lmim object s follows the
            % matlab convention (as for double array) for new
            % entries in array!
            p.lmim=s;
            h=p.lmim;
            h(end+1)=p.cost;

            var_old = p.var;
            bl_old = p.bl;

            %adjust varlist
            [p.var,p.bl]=vars_append(s, var_old, bl_old);            
        end
        function p = set.cost (p, s)
            p.cost = s;
            h=p.lmim;
            h(end+1)=p.cost;
            var_old = p.var;
            bl_old = p.bl;

            %adjust varlist
            [p.var,p.bl]=vars_append(s, var_old, bl_old);            

        end
        function p = set.fac (p, s)
            p.fac = s;            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %instantiaton can take any arguments
        %all other methods' first argument must be in class
        function so=lmis(p,s,c)
            % function so=lmis(p,s,c)
            %
            % If p is lmis and c double: Adjoin (s,c) to (p.lmim,p.fac)
            %              and c char:   Set p.cost to s.
            % If p is lmim and s double: Generate so with (so.lmim,so.fac)=(p,s)
            %                  s char:   Generate so with so.cost=s

            % generate new (pn,sn,cn): lmis pn, lmim sn, char or double cn
                if nargin==0;
                return
            end
            so.lmim=lmim.empty;
            so.fac=lmival.empty;
            sn=lmim;
            cn=[];
            if isa(p,'lmis')
                so=p;
                if nargin>=2;
                    sn=s;                    
                end
                if nargin>=3
                    cn=c;
                end
            else
                if ~isa(p,'lmim') && ~iscell(p)
                    error('First argument must be a single lmim.')
                else
                    % sn=lmim(p);
                    sn = p;
                    if nargin>2;
                        error('Only two arguments allowed.')
                    end
                    if nargin==2;
                        cn=s;
                    end
                end
            end
            %assign cost if cn is char
            if isa(cn,'char')
                [k,m]=size(sn.A);
                if k>1 | m>1
                    error('Cost map must be of dimension 1x1.')
                else
                    %cost map always replaced, no arrays allowed.
                    so.cost=sn;
                end
            elseif iscell(sn)
                lsn = length(sn);
                lf = length(so.fac);
                valcn = lmival(cn);

                for i = lsn:-1:1
                    so.lmim(lf + i) = sn{i};
                    so.fac(lf + i) = valcn;
                end
            else
                %otherwise adjoin sn to so.lmim and outer factor if
                %available
                so.lmim(end+1)=sn;
                so.fac(end+1)=lmival(cn);
                check(so)                
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
