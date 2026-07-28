classdef  opt_system_periodic_orbit < opt_system
    %OPT_SYSTEM_PERIODIC_ORBIT interconnection of network and operators
    

    
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
    % related together by orbit matrices M describing the symmetry
    %
    % Example:
    %  A(k)  = [Mx 0]^-k [A(0)  Bw(0)]  [Mx 0]^k
    %  Cz(k) = [0 Mz]    [Cz(k) Dzw(0)] [0 Mw]
    
    %make this extremely simple:
    %everything shuffles by M
    %and M is an orthogonal matrix

    properties
        M = [];  %symmetry matrix, should be orthogonal
        order = 0;  %M^order = M
    end
    
    methods
        function obj = opt_system_periodic_orbit(op, P, K, M, bind, tracking)
            %OPT_SYSTEM_PERIODIC constructor            
            if nargin < 6
                s = length(op);
                bind = 1:s;            
            end

            if nargin < 7
                 tracking = [];
            end

            obj@opt_system(op, P, K, bind, tracking)
            obj.M = M;
            obj.type = 'periodic_orbit';


            %order of the system
            orth_gap = norm(M' * M - eye(size(M)));
            if orth_gap > 1e-9 
                throw('Periodic Orbit: symmetry matrix M is not orthogonal')
            end


            Rc = M;
            order = 0;
            order_max = 1000;
            while (norm(Rc - eye(size(M))) > 1e-9) && (order < order_max)
                order = order + 1;
                Rc = Rc * M;
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

            if isnumeric(param) && isempty(param)
                param = struct('mode', 1);
            end
            Pcurr = obj.P.P;

            k = mod(param.mode-1, obj.order+1);           


            [n, m] = size(Pcurr.B);
            p = size(Pcurr.C, 1);
            c = size(obj.M, 1);

            Rxk = kron(eye(n/c), obj.M)^k;
            Ryk = kron(eye(p/c), obj.M)^k;
            Ruk = kron(eye(m/c), obj.M)^k;

            Pcurr.A = (Rxk) \ Pcurr.A * Rxk;
            Pcurr.B = (Rxk) \ Pcurr.B * Ruk;
            Pcurr.C = (Ryk) \ Pcurr.C * Rxk;
            Pcurr.D = (Ryk) \ Pcurr.D * Ruk;

        end

        function Kcurr = get_K(obj, param)
            %GET_K get the controller K
            if isa(obj.K, 'genplant')
                Kcurr = obj.K.ss;
            else
                Kcurr = obj.K;
            end

            if isnumeric(param) && isempty(param)
                param = struct('mode', 1);
            end


            if ~isempty(Kcurr)

                k = mod(param.mode-1, obj.order+1);          

                [n, m] = size(Kcurr.B);
                p = size(Kcurr.C, 1);
                c = size(obj.M, 1);

                Rxk = kron(eye(n/c), obj.M)^k;
                Ryk = kron(eye(m/c), obj.M)^k;
                Ruk = kron(eye(p/c), obj.M)^k;

                Kcurr.A = (Rxk) \ Kcurr.A * Rxk;
                Kcurr.B = (Rxk) \ Kcurr.B * Ruk;
                Kcurr.C = (Ryk) \ Kcurr.C * Rxk;
                Kcurr.D = (Ryk) \ Kcurr.D * Ruk;
            end
        end
 
        %% build the plant for the LMIs       
        function plant_rot = rotate_plant(obj, plant, direction)

            %rotate_plant: apply the periodic-orbit rotation to the
            %time-varying system, producing an LTI system
            if nargin < 3
                direction = 1;
            end

            c = size(obj.M, 1);
            n = size(plant.A, 1);
            Rkron = kron(eye(n/c), obj.M)^(direction);

            plant_rot = plant;
            if isa(plant, 'genplant')
                plant_rot.P.A = Rkron * plant.P.A;
                plant_rot.P.B = Rkron * plant.P.B;
            else
                plant_rot.A = Rkron * plant.A;
                plant_rot.B = Rkron * plant.B;
            end


        end

        function [alg_psi, iqc_op, alg_loop] = build_plant(obj, iqc_data, rho)
            %BUILD_PLANT: form the plant to be used for analysis
            %or synthesis
            %
            %Args:
            %   iqc_data: from manager.iqc_op_all, information about the
            %             operator iqc descriptions
            %   rho: exponential weighting
            %
            %Returns:
            %   alg_psi:    plant with filters (psi)
            %   iqc_op:     iqcs for the robust uncertainties
            %   alg_loop:   plant without filters, but after loop
            %               transformation (should be stable)

            if nargin < 3
                rho = 1;
            end



            %get the plant and the IQCs.

            if strcmp(iqc_data.task, 'analysis');
                alg0 = obj.get_alg([]);      
            else
                alg0 = obj.P.P;                
            end

            %rotate the plant
            if (iqc_data.rotate == 0)
                alg = alg0;
            else
                alg = obj.rotate_plant(alg0);
            end
            

            %repeat this call multiple times for switched systems. This
            %function will be overloaded, whereas build_plant_single will
            %stay the same.
            [alg_psi, iqc_op, alg_loop] = build_plant_single(obj, alg, iqc_data);



        end 


        %% exports
        function sys_per = export_periodic(obj)
            %EXPORT_PERIODIC export the periodic-orbit as a periodic system
            %explicitly list all subsystems
            %
            %Returns:
            %   sys_per (opt_system_periodic): a periodic system
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

            Pper_poly = genplant_poly(Pper);
            
            sys_per = opt_system_periodic(obj.op, Pper_poly, Kper, obj.bind, obj.tracking);
        end

        function sys_lift = periodic_lift(obj)
            %PERIODIC_LIFT form a periodic LTI lift of the system
            %create an equivalent LTI system
            %
            %Returns:
            %   sys_lift (opt_system): an LTI system

            sys_per = obj.export_periodic();

          sys_lift = sys_per.periodic_lift();

        end



        


    end
end

