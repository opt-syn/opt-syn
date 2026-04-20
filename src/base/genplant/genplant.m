classdef genplant
    %genplant a generalized plant
    
    
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
        function obj = genplant(P, n)
            %N Construct an instance of this class
            %   Detailed explanation goes here

            obj.P = P;
            if isnumeric(P)
                %static system
                [obj.ny, obj.nu] = size(P);
            else                     
                obj.s = n.s;
                obj.nz = n.nz;
                obj.nw = n.nw;
                obj.ny = n.ny;
                obj.nu = n.nu;            
    
                if isfield(n, 'zp')
                    obj.nzp = n.nzp;            
                end
                if isfield(n, 'zw')
                    obj.nwp = n.nwp;            
                end
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

        %% extract matrices       
        function B = Bw(obj)
            B0 = obj.B;
            B = B0(:, obj.index_w());
        end
        
        function B = Bu(obj)
            B0 = obj.B;
            B = B0(:, obj.index_u());
        end

        function C = Cz(obj)
            C0 = obj.C;
            C = C0(obj.index_z(), :);
        end
        
        function C = Cy(obj)
            C0 = obj.C;
            C = C0(obj.index_y(), :);
        end



        function D = Dzw(obj)
            %oracle output to oracle input
            iw = obj.index_w();
            iz = obj.index_z();

            D0 = obj.P.D;
            D = D0(iz, iw);
        end

        function D = Dzu(obj)
            %controller output to oracle input
            iu = obj.index_u();
            iz = obj.index_z();
            
            D0 = obj.P.D;
            D = D0(iz, iu);
        end

        function D = Dyw(obj)
            %oracle output to controller input
            iw = obj.index_w();
            iy = obj.index_y();

            D0 = obj.P.D;
            D = D0(iy, iw);
        end

        function D = Dyu(obj)
            %controller output to controller input
            iu = obj.index_u();
            iy = obj.index_y();

            D0 = obj.P.D;
            D = D0(iy, iu);
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

        function [Aa, B1, B2, C1, D11, D12, C2, D21, D22] = ss_zy_wu(obj)
            %get plant matrices for the [wu] -> [zy] subsytsem
            
            Aa = obj.A;
            B1 = obj.Bw;
            B2 = obj.Bu;
            C1 = obj.Cz;
            C2 = obj.Cy;
            D11 = obj.Dzw;
            D12 = obj.Dzu;
            D21 = obj.Dyw;
            D22 = obj.Dyu;

            
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

            P_diag = blkdiag(obj.P, b2.P);
            nin = obj.nw + obj.nwp + obj.nu;
            nout = obj.nz + obj.nzp + obj.ny;


            ind_in = [1:obj.nw, nin + (1:b2.nw), ...
                      obj.nw + (1:obj.nwp), nin + b2.nw + (1:b2.nwp), ...
                      obj.nw + obj.nwp + (1:obj.nu), nin + b2.nw + b2.nwp + (1:b2.nu)];

            ind_out = [1:obj.nz, nout + (1:b2.nz), ...
                      obj.nz + (1:obj.nzp), nout + b2.nz + (1:b2.nzp), ...
                      obj.nz + obj.nzp + (1:obj.ny), nout + b2.nz + b2.nzp + (1:b2.ny)];

            b_out.P = P_diag(ind_out, ind_in);
        end



        function b_out = lft(obj, b2)
            %LFT feedback interconnection of obj and plant bt
            %along common channels (u, y)

            b_out = obj;
            if isa(b2, 'genplant')
             
                
                b_out.nw = obj.nw;
                b_out.nwp = obj.nwp + b2.nwp;
                b_out.nz = obj.nz;
                b_out.nzp = obj.nzp + b2.nzp;
                b_out.nu = b2.nu;
                b_out.ny = b2.ny;
                b_out.s = obj.s + b2.s;
    
    
                b_out.P = lft(obj.P, b2.P, obj.nu, obj.ny);
            else
                [nu2, ny2] = size(b2.D);
                b_out.P = lft(obj.P, b2, nu2, ny2);
                b_out.nu = b_out.nu - nu2;
                b_out.ny = b_out.ny - ny2;
            end

            
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


        function n = dump_dim(obj)
            n = struct('nw', obj.nw, 'nwp', obj.nwp, ...
                'nu', obj.nw, 'ny', obj.ny, ...
                'nz', obj.nz, 'nzp', obj.nzp, ...
                's', obj.s);
        end
        %% performance inputs and outputs

        function [obj, iwp] = add_oracle_input(obj, ind_w, ind_z)

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

            iwp = obj.nwp + (1:(nwpnew + nzpnew));           
            obj.nwp = obj.nwp + nwpnew + nzpnew;           
        end
    
        function [obj, izp] = perf_output_w(obj, ind_w)
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

            izp = obj.nzp + (1:nnew);
            obj.nzp = obj.nzp + nnew;
            

        end

        function [obj, izp] = perf_output_opt(obj, c)
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

            izp = obj.nzp + (1:c);
            obj.nzp = obj.nzp + c;

        end

        function [obj, izp] = perf_output_z(obj, ind_z)
            %PERF_OUTPUT_Z: add performance to track the z output
            A = obj.P.A;
            B = obj.P.B;
            C = obj.P.C;
            D = obj.P.D;
            
            nnew = length(ind_z);

            Ctop = C([obj.index_z(), obj.index_zp()], :);
            Dtop = D([obj.index_z(), obj.index_zp()], :);
            Cbot = C([obj.index_y()], :);
            Dbot = D([obj.index_y()], :);

            n = length(A);
            
            
            Ez = full(sparse(1:nnew, ind_z, ones(nnew, 1), ...
                length(ind_z), nnew + obj.nzp + obj.ny));


            Czp = Ez * C;
            Dzp = Ez * D;
            

            Cnew = [Ctop; Czp; Cbot];
            Dnew = [Dtop; Dzp; Dbot];

            obj.P = ss(A, B, Cnew, Dnew, 1);

            izp = obj.nzp + (1:nnew);
            obj.nzp = obj.nzp + nnew;

        end
        

        function [obj, izp] = perf_output_con(obj, c, ind_z)
            
            %PERF_OUTPUT_CON: add performance to track the consensus output
            % norm(z)^2 (with z* = 0 by regulation)
            if nargin == 1
                c = 1;
            end
            if nargin == 2
                ind_z = 1:obj.nz;
            end

                        A = obj.P.A;
            B = obj.P.B;
            C = obj.P.C;
            D = obj.P.D;
            
            nnew = length(ind_z);

            Ctop = C([obj.index_z(), obj.index_zp()], :);
            Dtop = D([obj.index_z(), obj.index_zp()], :);
            Cbot = C([obj.index_y()], :);
            Dbot = D([obj.index_y()], :);

            n = length(A);
            
            %TODO: bug here.
            Ez = full(sparse(1:nnew, ind_z, ones(nnew, 1), ...
                length(ind_z), nnew + obj.nzp + obj.ny));
            
            Iz = eye(obj.nz);
            Jz = ones(obj.nz, obj.nz)/ (obj.nz/c);

            Resz = blkdiag((Iz - Jz), eye(obj.nzp + obj.ny));
            
            Czp = (Ez * Resz) *  C;
            Dzp = (Ez * Resz) * D;
            

            Cnew = [Ctop; Czp; Cbot];
            Dnew = [Dtop; Dzp; Dbot];

            obj.P = ss(A, B, Cnew, Dnew, 1);

            izp = obj.nzp + (1:nnew);
            obj.nzp = obj.nzp + nnew;            
        end
    
    end
end

