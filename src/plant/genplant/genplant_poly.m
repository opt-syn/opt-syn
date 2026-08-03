classdef genplant_poly <  genplant & matlab.mixin.indexing.RedefinesBrace
    %GENPLANT_POLY a generalized plant defined over corners of a polytope.
    %Used for switched systems

    
    % properties
    %     %plant matrices
    % end

    methods (Access=public)
        function obj = genplant_poly(P_orig, n)
            %Constructor
            %
            %Args:
            %   P (cell of sdpss): state space systems for each subsystem
            %   n (struct): partition of the channels into [z, zp, y], [w, wp, u]

           
     
            
            if isempty(P_orig)
                n = [];
            else
                if iscell(P_orig)
                    P = P_orig;
                else
                    Nss = length(P_orig);
                    P = cell(Nss, 1);
                    for i = 1:Nss
                        P{i} = P_orig{i};
                    end            
                end
                n = P{1}.dump_dim();                
            end

            obj@genplant(P, n);


            
        end

        
        function obj = rhotrafo(obj, rho)
            % RHOTRAFO Apply an exponential discount (rho-transformation).
            %
            % Scales ``A`` and ``B`` by :math:`\rho^{-1}`, which corresponds
            % to an exponential weighting of the signals in discrete time.
            %
            % :param rho: Discount factor.
            % :type rho: double
            % :returns: The transformed plant (modified in place).
            % :rtype: genplant
            for i = 1:obj.Nss
                obj.P{i} = obj.P{i}.rhotrafo(rho);
            end
        end

        function b_out = lft(obj, b2)
            %LFT linear fractional transformation:
            %feedback interconnection of obj and plant b2
            %along common channels (u, y) in each subsystem 
            %obj star b2
            %
            %Args:
            %   b2 (genplant): the other plant
            %Returns
            %   b_out (genplant): the lft plant

            % Initialize the output for the linear fractional transformation
            b_out = obj; 
            for i = 1:obj.Nss

                if iscell(b2)
                    b_out.P{i} = lft(obj.P{i}, b2{i});
                elseif isa(b2, 'genplant_poly')
                    b_out.P{i} = lft(obj.P{i}, b2.P{i});
                else
                    b_out.P{i} = lft(obj.P{i}, b2);
                end                
            end

            b_out.nw = b_out.P{1}.nw;
            b_out.nwp = b_out.P{1}.nwp;
            b_out.nz = b_out.P{1}.nz;
            b_out.nzp = b_out.P{1}.nzp;
            b_out.nu = b_out.P{1}.nu;
            b_out.ny = b_out.P{1}.ny;
            b_out.s = b_out.P{1}.s;

        end

        function D = Dyu(obj)
            %controller output to controller input, direct feedthrough
            D = cellfun(@Dyu, obj.P);            
        end

        function Ao = A(obj)
            %cell of A matrices in state space systems
            Ao = cellfun(@A, obj.P, 'UniformOutput',false);            
        end
        function Bo = B(obj)
            %cell of B matrices in state space systems
            Bo = cellfun(@B, obj.P, 'UniformOutput',false);            
        end
        function Co = C(obj)
            %cell of C matrices in state space systems
            Co = cellfun(@C, obj.P, 'UniformOutput',false);            
        end
        function Do = D(obj)
            %cell of D matrices in state space systems
            Do = cellfun(@D, obj.P, 'UniformOutput',false);            
        end       

        function nxo = nx(obj)
            %number of states
            nxo = length(obj.P{1}.A);
        end

        function nss = Nss(obj)
            %number of subsystems
            nss = length(obj.P);
        end

        function P_out = ss(obj)
            %extract the state-space expression
            P_out = cellfun(@ss, obj.P);
        end

        function [Aa, B1, B2, C1, D11, D12, C2, D21, D22] = ss_zy_wu(obj, ind)
            %get plant matrices for the [wu] -> [zy] subsytsem
            %
            %Args:
            %   ind: the index of the subsystem to get information from
            
            [Aa, B1, B2, C1, D11, D12, C2, D21, D22] = obj.P{ind}.ss_zy_wu();
        end
        %% overloads

        function b_out = blkdiag(obj, b2)
            %block-diagonal of two plants, comporting with the indexing scheme
            %
            %Args:
            %   b2 (genplant): the other plant
            %Returns
            %   b_out (genplant): the block diagonal plant

            %interleave the indices properly

            b_out = obj;
            b_out.nw = obj.nw + b2.nw;
            b_out.nwp = obj.nwp + b2.nwp;
            b_out.nz = obj.nz + b2.nz;
            b_out.nzp = obj.nzp + b2.nzp;
            b_out.nu = obj.nu + b2.nu;
            b_out.ny = obj.ny + b2.ny;
            b_out.s = obj.s + b2.s;

            for i = 1:obj.Nss
                b_out.P{i} = blkdiag(obj.P{i}, b2.P{i});
            end
        end

        
        %% performance inputs and outputs
        function obj = add_oracle_input(obj, ind_w, ind_z)
            %add external inputs at the operator F
            %
            %:math:`w + \delta w \in F(z + dz)`
            %Args:
            %   ind_w: at the input of the operator
            %   ind_z: at the output of the operator
            %Return:
            %   iwp: new performance input indices 

            %
            %Does not add extra outputs
            
            nwpnew = length(ind_w);
            nzpnew = length(ind_z);

            obj.P = cellfun(@(p) p.add_oracle_input(ind_w, ind_z), obj.P, 'UniformOutput', false);

            obj.nwp = obj.nwp + nwpnew + nzpnew;           
        end
    
        function obj = perf_output_w(obj, ind_w)
            %PERF_OUTPUT_W: add performance to track the w output
            %
            %Args:
            %   ind_w: indices of the input of the operator
            %Return:
            %   izp: new performance output indices

            nnew = length(ind_w);
            obj.P = cellfun(@(p) p.perf_output_w(ind_w), obj.P);
            obj.nzp = obj.nzp + nnew;

        end

        function obj = perf_output_opt(obj, c)
            %PERF_OUTPUT_OPT add performance to track the w output
            %condition,   
            %Args:
            %   ind_w: indices of the input of the operator
            %Return:
            %   izp: new performance output indices
            %

            % :math:`z_p = \sum_{i=1}^s w^i`
            if nargin == 1
                c = 1;
            end            
            obj.P = cellfun(@(p) p.perf_output_opt(c), obj.P, 'UniformOutput', false);
            
            obj.nzp = obj.nzp + c;

        end

        function obj = perf_output_z(obj, ind_z)
            %PERF_OUTPUT_Z: add performance to track the z output
            %
            %Args:
            %   ind_z: indices of the output of the operator
            %Return:
            %   izp: new performance output indices

            nnew = length(ind_z);

            obj.P = cellfun(@(p) p.perf_output_z(ind_z), obj.P, 'UniformOutput', false);

            obj.nzp = obj.nzp + nnew;

        end
        

        function obj = perf_output_con(obj, c, ind_z)
            %PERF_OUTPUT_CON: add performance to track the consensus output
            % `z_p^i = z^i - \text{average}(z)`
            %
            %Args:
            %   c: dimension of the coordinate/kronecker lift
            %   iz: indices of input of operators
            %Return:
            %   izp: new performance output indices

            if nargin == 1
                c = 1;
            end
            if nargin == 2
                ind_z = 1:obj.nz;
            end


            
            nnew = length(ind_z);

            obj.P = cellfun(@(p) p.perf_output_con(c, ind_z), obj.P, 'UniformOutput', false);

            obj.nzp = obj.nzp + nnew;            
        end

        function [obj, iwp, izp] = perf_ergodic(obj, Nw)
            %PERF_ERGODIC inputs and outputs for ergodic convergence
            %
            %Args:
            %   Nw: given consensus matrix
            %Return:
            %   iwp: new performance input indices
            %   izp: new performance output indices

            for i = 1:obj.Nss
                [obj.P{i}, iwp, izp] = obj.P{i}.perf_ergodic(Nw);
            end
            obj.nwp = obj.nwp + length(iwp);
            obj.nzp = obj.nzp + length(izp);
        end

    end

    methods (Access=protected)
        
        
        %https://de.mathworks.com/help/matlab/ref/matlab.mixin.indexing.redefinesbrace-class.html
        function P_out = braceReference(obj,ind)
           P_out = obj.P.(ind);
        end

        function P_out = braceListLength(obj,indexOp,indexContext)
           P_out = listLength(obj.P,indexOp,indexContext);
        end

        function obj = braceAssign(obj,indexOp,varargin)
            if isscalar(indexOp)
                [obj.P.(indexOp)] = varargin{:};
                return;
            end
            [obj.P(indexOp)] = varargin{:};
        end
    end

end

