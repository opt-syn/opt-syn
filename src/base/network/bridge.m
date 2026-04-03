classdef bridge
    %BRIDGE a network sitting between the oracle F and the controller K
    
    
    properties
        %plant matrices
        P;
        
        %indexing [iz, izp, iy], [iw, iwp, iu]
        s = 0;
        nz = 0; %input to operators (from network)        
        nzp =0; %input to performance channel (from network)
        ny =0;  %input to controller (from network)


        nw =0;  %output of operators (to network)
        nwp =0; %output of performance channel (to network)
        nu =0;  %output of controller (to network)        
           
    end
    
    methods
        function obj = bridge(P, n)
            %N Construct an instance of this class
            %   Detailed explanation goes here
            obj.P = P;
%             [obj.A, obj.B, obj.C, obj.D] = ssdata(P);            
            obj.s = n.s;
            obj.nz = n.nz;
            obj.nw = n.nw;
            obj.ny = n.ny;
            obj.nu = n.nu;            

            if isfield(n, 'zp')
                obj.nzp = n.nzp;
            else
                obj.nzp = 0;
            end
            if isfield(n, 'zw')
                obj.nwp = n.nwp;
            else
                obj.nwp = 0;
            end
            
        end

        %indexers
        function u_ind = index_u(obj)
            u_ind = obj.nw + obj.nwp + (1:obj.nu);
        end

        function wp_ind = index_wp(obj)
            wp_ind = obj.nw + (1:obj.nwp);
        end

        function w_ind = index_w(obj)
            w_ind = 1:obj.nw;
        end

        function y_ind = index_y(obj)
            y_ind = obj.nz + obj.nzp + (1:obj.ny);
        end

        function zp_ind = index_zp(obj)
            zp_ind = obj.nz + (1:obj.nzp);
        end

        function z_ind = index_z(obj)
            z_ind = 1:obj.nz;
        end

        function wr_ind = index_notu(obj)
            %TODO: expand with more inputs
            wr_ind = 1:(obj.nw + obj.nwp);
        end

        %extract matrices       
        function D = Dyu(obj)
            %get the direct feedthrough matrix
            iu = obj.index_u();
            iy = obj.index_y();

            D = obj.D(iy, iu);
        end

        function Ao = A(obj)
            Ao = obj.P.A;
        end
        function Bo = B(obj)
            Bo = obj.P.B;
        end
        function Co = C(obj)
            Co = obj.P.C;
        end
        function Do = D(obj)
            Do = obj.P.D;
        end

        function nxo = nx(obj)
            nxo = length(obj.P.A);
        end

        function P_out = ss(obj)
            %extract the state-space expression
            P_out = obj.P;
        end

        function P_out = tf(obj)
            P_out= ss2tf(obj.ss());
        end

        function obj = lift(obj, d)
            %lift by a kronecker operation with the identity            
            
            Ad = kron(obj.P.A, eye(d));
            Bd = kron(obj.P.B, eye(d));
            Cd = kron(obj.P.C, eye(d));
            Dd = kron(obj.P.D, eye(d));

            obj.P = ss(Ad, Bd, Cd, Dd, 1);

            obj.nz = obj.nz * d;
            obj.nzp = obj.nzp * d;
            obj.nw = obj.nw * d;
            obj.nwp = obj.nwp * d;
            obj.nu = obj.nu * d;
            obj.ny = obj.ny * d;
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

            B = obj.B;
            D = obj.D;
            if nwpnew 
                Ew = full(sparse(ind_w, 1:nwpnew, ones(nwpnew, 1), obj.nw, nwpnew));
                Bwnew = B(:, obj.index_w) * Ew;
                Dwnew = D(:, obj.index_w) * Ew;
%                 Bwnew = 0 * [Ew; zeros(nzpnew, nwpnew)];
%                 Dwnew = 1 * [Ew; zeros(nzpnew, nwpnew)];
            else
                Bwnew = [];
                Dwnew = [];
            end
            if nzpnew
                Ez = full(sparse(ind_z, 1:nzpnew, ones(nzpnew, 1), obj.nz, nzpnew));
                Bznew = [];
                Dznew = 1 * [ Ez; zeros(nzpnew, obj.nz)];
            else
                Bznew = [];
                Dznew = [];
            end

            
            

            B_left = B(:, [obj.index_w, obj.index_wp]);
            B_right = B(:, [obj.index_u]);
            D_left = D(:, [obj.index_w, obj.index_wp]);
            D_right = D(:, [obj.index_u]);

            Aold = obj.P.A;            
            Bnew = [B_left, Bwnew, Bznew, B_right];
            Cold = obj.P.C;
            Dnew = [D_left, Dwnew, Dznew, D_right];

            obj.P = ss(Aold, Bnew, Cold, Dnew, 1);

            obj.nwp = obj.nwp + nwpnew + nzpnew;           
        end
    
        function obj = perf_output_w(obj, ind_w)
            %PERF_OUTPUT_W: add performance to track the w output
            A = obj.P.A;
            B = obj.P.B;
            C = obj.P.C;
            D = obj.P.D;
            
            nnew = length(ind_w);

            Ctop = C([obj.index_z(), obj.index_zp()], :);
            Dtop = D([obj.index_z(), obj.index_zp()], :);
            Cbot = C([obj.index_y()], :);
            Dbot = D([obj.index_y()], :);

            n = length(A);
            Czp = zeros(nnew, n);
            Ez = full(sparse(ind_w, 1:nnew, ones(nnew, 1), n, nnew));

            Dzp = Ez;

            Cnew = [Ctop; Czp; Cbot];
            Dnew = [Dtop; Dzp; Dbot];

            obj.P = ss(A, B, Cnew, Dnew, 1);

            obj.nzp = obj.nzp + nnew;

        end

        function obj = perf_output_opt(obj, c)
            %PERF_OUTPUT_WSUM: add performance to track the optimality
            %condition: sum(1'w) = 0
            %
            if nargin == 1
                c = 1;
            end
            A = obj.P.A;
            B = obj.P.B;
            C = obj.P.C;
            D = obj.P.D;

            ind_w = obj.index_w();
            
            nnew = length(ind_w);

            Ctop = C([obj.index_z(), obj.index_zp()], :);
            Dtop = D([obj.index_z(), obj.index_zp()], :);
            Cbot = C([obj.index_y()], :);
            Dbot = D([obj.index_y()], :);

            n = length(A);
            s = obj.nw/c;
            Czp = zeros(c, n);


            Ezp = kron(ones(1, s), eye(c));
            Dzp = [Ezp, zeros(size(Ezp, 1), obj.nwp + obj.nu)];

            Cnew = [Ctop; Czp; Cbot];
            Dnew = [Dtop; Dzp; Dbot];

            obj.P = ss(A, B, Cnew, Dnew, 1);

            obj.nzp = obj.nzp + c;

        end

        function obj = perf_output_z(obj, ind_z)
            %PERF_OUTPUT_Z: add performance to track the z output
            A = obj.P.A;
            B = obj.P.B;
            C = obj.P.C;
            D = obj.P.D;
            
            nnew = length(ind_z);

            Ctop = C([obj.index_z(), obj.index_zp()], :);
            Dtop = C([obj.index_z(), obj.index_zp()], :);
            Cbot = C([obj.index_y()], :);
            Dbot = D([obj.index_y()], :);

            n = length(A);
            
            %TODO: bug here.
            Ez = full(sparse(ind_z, 1:nnew, ones(nnew, 1), length(ind_z), ...
                nnew + obj.nw));


            Czp = Ez * C;
            Dzp = Ez * D;
            

            Cnew = [Ctop; Czp; Cbot];
            Dnew = [Dtop; Dzp; Dbot];

            obj.P = ss(A, B, Cnew, Dnew, 1);

            obj.nzp = obj.nzp + nnew;

        end
        

        function obj = perf_output_con(obj)
            %PERF_OUTPUT_CON: add performance to track the consensus output
            % norm(z)^2 (with z* = 0 by regulation)
            obj = obj.perf_output_z(obj.index_z());
        end
    
    end
end

