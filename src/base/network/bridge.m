classdef bridge
    %N Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %plant matrices
        A;      
        B;          
        C;
        D;
        
        %indexing
        nz = 0;
        nzp =0;
        ny =0;
        nw =0;
        nwp =0;
        nu =0;
        nx = 0;
%         iz;     %input to operators (from network)        
%         izp;    %input to performance channel (from network)
%         iy;     %input to controller (from network)
% 
% 
%         iw;     %output of operators (to network)
%         iwp;    %output of performance channel (to network)
%         iu;     %output of controller (to network)        
    end
    
    methods
        function obj = bridge(P, n)
            %N Construct an instance of this class
            %   Detailed explanation goes here
            [obj.A, obj.B, obj.C, obj.D] = ssdata(P);            

            obj.nz = n.nz;
            obj.nw = n.nw;
            obj.ny = n.ny;
            obj.nu = n.nu;
            obj.nx = length(obj.A);

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

        function P_out = ss(obj)
            %extract the state-space expression
            P_out = ss(obj.A, obj.B, obj.C, obj.D, 1);
        end

        function P_out = tf(obj)
            P_out= ss2tf(obj.ss());
        end

        function obj = lift(obj, d)
            %lift by a kronecker operation with the identity            
            obj.P = kron(obj.P, eye(d));

            obj.nz = obj.nz * d;
            obj.nzp = obj.nzp * d;
            obj.nw = obj.nw * d;
            obj.nwp = obj.nwp * d;
            obj.nu = obj.nu * d;
            obj.ny = obj.ny * d;
        end
    end
end

