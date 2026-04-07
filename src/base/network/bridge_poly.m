classdef bridge_poly <  bridge & matlab.mixin.indexing.RedefinesBrace
    %BRIDGE_POLY a network sitting between the oracle F and the controller K

    
    % properties
    %     %plant matrices
    % end

    methods (Access=public)
        function obj = bridge_poly(P, n)
            %N Construct an instance of this class
            %   Detailed explanation goes here
            
           
            obj@bridge([], n);

            Nss = length(P);
            obj.P = cell(Nss, 1);
            for i = 1:Nss
                obj.P{i} = P{i};
            end
            
        end
        function D = Dyu(obj)
            %get the direct feedthrough matrix
            D = cellfun(@Dyu, obj.P);            
        end

        function Ao = A(obj)
            Ao = cellfun(@A, obj.P, 'UniformOutput',false);            
        end
        function Bo = B(obj)
            Bo = cellfun(@B, obj.P, 'UniformOutput',false);            
        end
        function Co = C(obj)
            Co = cellfun(@C, obj.P, 'UniformOutput',false);            
        end
        function Do = D(obj)
            Do = cellfun(@D, obj.P, 'UniformOutput',false);            
        end       

        function nxo = nx(obj)
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

        function P_out = tf(obj)
            P_out= ss2tf(obj.ss());
        end

        %% overloads

        function b_out = blkdiag(obj, b2)
            %block-diagonal of two bridges
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

            %ADD_ORACLE_INPUT: add external inputs at the oracle F
            %
            %z + dz \in F(w + dz) + dz
            %ind_w: at the input of the oracle
            %ind_z: at the output of the oracle

            %
            %Does not add extra outputs
            
            nwpnew = length(ind_w);
            nzpnew = length(ind_z);

            obj.P = cellfun(@(p) p.add_oracle_input(ind_w, ind_z), obj.P);

            obj.nwp = obj.nwp + nwpnew + nzpnew;           
        end
    
        function obj = perf_output_w(obj, ind_w)
            %PERF_OUTPUT_W: add performance to track the w output
            
            nnew = length(ind_w);
            obj.P = cellfun(@(p) p.perf_output_w(ind_w), obj.P);
            obj.nzp = obj.nzp + nnew;

        end

        function obj = perf_output_opt(obj, c)
            %PERF_OUTPUT_WSUM: add performance to track the optimality
            %condition: sum(1'w) = 0
            %
            if nargin == 1
                c = 1;
            end            
            obj.P = cellfun(@(p) p.perf_output_opt(c), obj.P);
            
            obj.nzp = obj.nzp + c;

        end

        function obj = perf_output_z(obj, ind_z)
            %PERF_OUTPUT_Z: add performance to track the z output
            
            nnew = length(ind_z);

            obj.P = cellfun(@(p) p.perf_output_z(ind_z), obj.P);

            obj.nzp = obj.nzp + nnew;

        end
        

        function obj = perf_output_con(obj, c, ind_z)
            
            %PERF_OUTPUT_CON: add performance to track the consensus output
            % norm(z)^2 (with z* = 0 by regulation)
            if nargin == 1
                c = 1;
            end
            if nargin == 2
                ind_z = 1:obj.nz;
            end


            
            nnew = length(ind_z);

            obj.P = cellfun(@(p) p.perf_output_con(c, ind_z), obj.P);

            obj.nzp = obj.nzp + nnew;            
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

