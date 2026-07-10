classdef  opt_system_periodic < opt_system_switched
    %OPT_SYSTEM_PERIODIC interconnection of network and operators
    %
    %a periodic system: repeated and predictable cycle evaluation   
    %
    % w(k) \in F(z(k))
    %
    % [x(k+1)] = [A(k)    Bw(k)    Bwp(k)   Bu(k)  ][x(k)]   state transition
    % [z(k)  ] = [Cz(k)   Dzw(k)   Dzwp(k)  Dzu(k) ][w(k)]   output to oracle
    % [zp(k) ] = [Czp(k)  Dzpw(k)  Dzpwp(k) Dzpu(k)][wp(k)]  output to performance
    % [y(k)  ] = [Cy(k)   Dyw(k)   Dywp(k)  Dyu(k) ][u(k)]   output to controller
    %
    %All matrices periodic: A(k) = A(k+T) for some known time T
    %(T = number of subsystems)

    
    methods
        function obj = opt_system_periodic(op, P, K, bind, tracking)
            %OPT_SYSTEM_PERIODIC constructor            
            if nargin < 4
                s = length(op);
                bind = 1:s;            
            end

            if nargin < 5
                 tracking = [];
            end

            Nss = max(P.Nss);
            adj = circshift(eye(Nss), -1);
            obj@opt_system_switched(op, P, K, adj, bind, tracking)
            obj.type = 'periodic';
        end        


        function [Sbeta, Rbeta] = get_tracked_opt(obj)
            %GET_TRACKED_OPT get the tracked position of the optimal
            %solution. allow for time-varying exosystems (periodic),
            %represented by a cell

            [Sbeta0, Rbeta0] = get_tracked_opt@opt_system_switched(obj);

            if iscell(Sbeta0)
                Sbeta = Sbeta0;
                Rbeta = Rbeta0;
            else
                Sbeta = cell(obj.P.Nss, 1);
                Rbeta = cell(obj.P.Nss, 1);

                for i = 1:obj.P.Nss
                    Sbeta{i} = Sbeta0;
                    Rbeta{i} = Rbeta0;
                end
            end
 
        end

        function sys_lift = periodic_lift(obj)

            %PERIODIC_LIFT form a periodic LTI lift of the system
            %create an equivalent LTI system

            sys_per = periodic_lift(obj.P.P);
            


            %dimension counters
            Nss = obj.Nss;
            n0= obj.P.dump_dim();

            %all dimensions
            nin = n0.nw + n0.nwp + n0.nu;
            nout = n0.nz + n0.nzp + n0.ny;

            %index the inputs
            i0_in = 1:(Nss*nin);
            ind_in = reshape(1:(Nss*nin), nin, []);
            ind_w = reshape(ind_in(1:n0.nw, :), [], 1);
            ind_wp = reshape(ind_in(n0.nw + (1:n0.nwp), :), [], 1);
            ind_u = reshape(ind_in(n0.nw + n0.nwp + (1:n0.nu), :), [], 1);

            % perm_in = i0_in([ind_w; ind_wp; ind_u]);

            %index the outputs
            i0_out = 1:(Nss*nout);
            ind_out = reshape(1:(Nss*nout), nout, []);
            ind_z = reshape(ind_out(1:n0.nz, :), [], 1);
            ind_zp = reshape(ind_out(n0.nz + (1:n0.nzp), :), [], 1);
            ind_y = reshape(ind_out(n0.nz + n0.nzp + (1:n0.ny), :), [], 1);

            % perm_out = i0_out([ind_z; ind_zp; ind_y]);

            % sys_perm = sys_per(perm_out, perm_in);

            sys_perm = sys_per([ind_z; ind_zp; ind_y], [ind_w; ind_wp; ind_u]);

            %track the dimensions
            n = n0;
            %outputs
            n.nz = n.nz * Nss;
            n.nzp = n.nzp * Nss;
            n.ny = n.ny * Nss;

            %inputs
            n.nw = n.nw * Nss;
            n.nwp = n.nwp * Nss;
            n.nu= n.nu * Nss;
            

            %lift of the controller

            if all(cellfun(@isempty, obj.K))
                sys_K = [];
            else
                sys_K = periodic_lift(obj.K);
            end

            

            P_lift = genplant(sys_perm, n);
            %lift of the trackers

            if ~isempty(obj.tracking)

            end

            bind_lift = repmat(obj.bind, 1, Nss);
            sys_lift = opt_system(obj.op, P_lift, sys_K, bind_lift, obj.tracking);

            % sys_lift = genplant(sys_perm, n);
        end

        function tp = get_type(obj)
            %get the type of the switched system
            %is periodic!
            tp = obj.type;
        end

        function mode_next = next_mode(obj, mode)
            %next mode in the switching sequence
            nss = obj.P.Nss;
            mode_next = 1+ mod(mode, nss);
        end
    end
end

