classdef  opt_system_periodic_orbit < opt_system
    %OPT_SYSTEM_PERIODIC_ORBIT interconnection of network and operators
    %
    %a periodic system: repeated and predictable cycle evaluation  
    %
    % orbit: the periodicity is highly structured in a symmetric manner
    %
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
    %
    % related together by orbit matrices R describing the symmetry
    %
    % Example:
    %  A(k)  = [Rx 0]^-k [A(0)  Bw(0)]  [Rx 0]^k
    %  Cz(k) = [0 Rz]    [Cz(k) Dzw(0)] [0 Rw]
    
    %make this extremely simple:
    %everything shuffles by R
    %and R is an orthogonal matrix

    properties
        R = [];
        order = 0;
    end
    
    methods
        function obj = opt_system_periodic_orbit(op, P, K, R, bind, tracking)
            %OPT_SYSTEM_PERIODIC constructor            
            if nargin < 6
                s = length(op);
                bind = 1:s;            
            end

            if nargin < 7
                 tracking = [];
            end

            obj@opt_system(op, P, K, bind, tracking)
            obj.R = R;
            obj.type = 'periodic_orbit';


            %order of the system
            orth_gap = norm(R' * R - eye(size(R)));
            if orth_gap > 1e-9 
                throw('Periodic Orbit: symmetry matrix R is not orthogonal')
            end


            Rc = R;
            order = 0;
            order_max = 1000;
            while (norm(Rc - eye(size(R))) > 1e-9) && (order < order_max)
                order = order + 1;
                Rc = Rc * R;
            end

            obj.order = order;

            
        end 

        function ns = Nss(obj)
            ns = obj.order+1;
        end


        function tp = get_type(obj)
            %get the type of the switched system
            %is periodic!
            tp = obj.type;
        end

        function mode_next = next_mode(obj, mode)
            %next mode in the switching sequence
            nss = obj.Nss;
            mode_next = 1+ mod(mode, nss);
        end

        %fetch attributes
        function Pcurr = get_P(obj, param)
            %GET_P get the network P
            Pcurr = obj.P.P;

            k = mod(param.mode-1, obj.order+1);           


            [n, m] = size(Pcurr.B);
            p = size(Pcurr.C, 1);
            c = size(obj.R, 1);

            Rxk = kron(eye(n/c), obj.R)^k;
            Ryk = kron(eye(p/c), obj.R)^k;
            Ruk = kron(eye(m/c), obj.R)^k;

            Pcurr.A = (Rxk) \ Pcurr.A * Rxk;
            Pcurr.B = (Rxk) \ Pcurr.B * Ruk;
            Pcurr.C = (Ryk) \ Pcurr.C * Rxk;
            Pcurr.D = (Ryk) \ Pcurr.D * Ruk;

        end

        function [Aa, B1, B2, C1, D11, D12, C2, D21, D22] = ss_zy_wu(obj, param)
            %get state space matrices at the current parameter values
            

            [Aa, B1, B2, C1, D11, D12, C2, D21, D22] = obj.P.ss_zy_wu();

            k = mod(param.mode-1, obj.order+1);          


            [nx, nw] = size(B1);
            [ny, nu] = size(C2);
            [nz] = size(D11, 1);
            

            Rxk = kron(eye(nx/c), obj.R)^k;
            Ryk = kron(eye(ny/c), obj.R)^k;
            Ruk = kron(eye(nu/c), obj.R)^k;
            Rzk = kron(eye(nx/c), obj.R)^k;
            Rwk = kron(eye(ny/c), obj.R)^k;   


            %perform the operation
            Aa = Rxk \ Aa * Rxk;
            B1 = Rxk \ B1 * Rwk;
            B2 = Rxk \ B2 * Ruk;

            C1 = Rzk \ C1 * Rxk;
            D11 = Rzk \ D11 * Rwk;
            D12 = Rzk \ D12 * Ruk;

            C2 = Ryk \ C2 * Rxk;
            D21 = Ryk \ D11 * Rwk;
            D22 = Ryk \ D21 * Ruk;



        end

        function [Sbeta, Rbeta] = get_tracked_opt(obj, param)
            %GET_TRACKED_OPT get the tracked position of the optimal
            %solution
            c = size(obj.R);
            k = mod(param.mode-1, obj.order+1);          

            
            if isempty(obj.tracking)
                Sbeta = eye(c);
                Rbeta = eye(c);
            else
                [ne, nd] = size(obj.tracking.Rbeta);

                Rdk = kron(eye(nd/c), obj.R)^k;
                Rek = kron(eye(ne/c), obj.R)^k;

                Sbeta = Rdk \ obj.tracking.Sbeta * Rdk;
                Rbeta = Rek \ obj.tracking.Rbeta * Rdk;
            end
        end

        function sys_per = export_periodic(obj)
            %EXPORT_PERIODIC export the periodic-orbit as a periodic system
            %explicitly list all subsystems
            n = obj.P.dump_dim();

            nss = obj.Nss;
            Pper = cell(1, nss);
            Kper = cell(1, nss);

            for i = 1:nss
                param_curr = struct('mode', i);
                Pper_sys = obj.get_P(param_curr);
                Kper_sys = obj.get_K(param_curr);


                Pper{i} = genplant(Pper_sys, n);
                Kper{i} = Kper_sys;
            end

            Pper_poly = genplant_poly(Pper, n);
            
            sys_per = opt_system_periodic(obj.op, Pper_poly, Kper, obj.bind, obj.tracking);
        end

        function sys_lift = periodic_lift(obj)
            %PERIODIC_LIFT form a periodic LTI lift of the system
            %create an equivalent LTI system

            sys_per = obj.export_periodic();

            sys_lift = sys_per.periodic_lift();

        end



        function Kcurr = get_K(obj, param)
            %GET_K get the controller K
            if isa(obj.K, 'genplant')
                Kcurr = obj.K.ss;
            else
                Kcurr = obj.K;
            end


            if ~isempty(Kcurr)

                k = mod(param.mode-1, obj.order+1);          
    
                [n, m] = size(Kcurr.B);
                p = size(Kcurr.C, 1);
                c = size(obj.R, 1);
    
                Rxk = kron(eye(n/c), obj.R)^k;
                Ryk = kron(eye(m/c), obj.R)^k;
                Ruk = kron(eye(p/c), obj.R)^k;
    
                Kcurr.A = (Rxk) \ Kcurr.A * Rxk;
                Kcurr.B = (Rxk) \ Kcurr.B * Ruk;
                Kcurr.C = (Ryk) \ Kcurr.C * Rxk;
                Kcurr.D = (Ryk) \ Kcurr.D * Ruk;
            end
        end


    end
end

