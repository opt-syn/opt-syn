classdef bridge
    %N Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        P;      %plant
        
        %indexing
        nz = 0;
        nzp =0;
        ny =0;
        nw =0;
        nwp =0;
        nu =0;
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
            obj.P = P;            

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

        %lft operations        

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

        function obj = lft_upper(obj, upper_part, ind)
            %close the loop as Q star P (from above) at the indices ind
            
        end

        function obj = lft_lower(obj, lower_part, ind_u, ind_y)
            %close the loop as  P star K (from below) at the indices ind
% 
%             if nargin < 3
%                 ind_u = 
%             end

        end
    end
end

