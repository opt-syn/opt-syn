classdef  opt_system_lpv_poly < opt_system_interface
    %OPT_SYSTEM_LPV_POLY interconnection of network and operators
    %
    %an lpv-polytopic system: gain-scheduled controller

    
    properties
        %polytopes: convex hulls
        par = [];  %polytope of parameter values 
        pair = []; %polytope of transitions [parameter, next_paramater] 

        % interp = []; %a vertex interpolator: returns a set of weights
        

        %vertex-scheduled controller
        %P: number of parameters
        %K: number of corners of the parameter polytope
       
    end

    methods
        function obj = opt_system_lpv_poly(op, P, K, pair, bind, tracking)
            %OPT_SYSTEM_PERIODIC constructor    


             if nargin < 5
                s = length(op);
                bind = 1:s;            
            end

            if nargin < 6
                 tracking = [];
            end


            obj@opt_system_interface(op, P, K, bind, tracking);
            %process the polytope set

            Npar = P.Nss;
            Ncorner = size(pair, 1);
            if size(pair, 2) == Npar
                %all to all transitions
                %duplicate the polytope values
                obj.par = pair;
                obj.pair = [kron(ones(Npar, 1), obj.par), kron(obj.par, ones(Npar, 1))];

                Npair = Npar^2;
            else
                obj.pair = pair;
                par_iso = obj.pair(:, 1:Npar);
                obj.par = unique(par_iso', 'rows')';
            end

            %also need a vertex-scheduled internal model (?)

            %vertex scheduled controller (default to the same controller
            %everywhere)
            if ~iscell(K)
                K0 = K;
                K = cell(Ncorner, 1);
                for i = 1:Ncorner
                    K{i} = K0;
                end
            end

            if ~isempty(obj.tracking) && ~iscell(obj.tracking.Sbeta)
                Sbeta0 = obj.tracking.Sbeta;
                Rbeta0 = obj.tracking.Rbeta;

                obj.tracking.Sbeta = cell(Npar, 1);
                obj.tracking.Rbeta = cell(Npar, 1);

                for i = 1:Npar
                    obj.tracking.Sbeta{i} = Sbeta0;
                    obj.tracking.Rbeta{i} = Rbeta0;
                end
            end
            % obj.interp = get_vertex_interp(obj.par);
            obj.type = 'lpv_poly';
        end        

        %% indexing the polytope
        function nss = Nss(obj)
            %number of parameters (subsystems)
            nss= obj.P.Nss;
        end

        function nss = Ncorner(obj)
            %number of parameters (subsystems)
            nss= size(obj.par, 1);
        end

        function nss = Npair(obj)
            %number of parameters (subsystems)
            nss= size(obj.pair, 1);
        end

        function dimn = nxn(obj)
            %nxn: number of states in network
            dimn = obj.P.nx;
        end

        function dimn = nzp(obj)
            %nzp: number of performance outputs
            dimn = obj.P.nzp;
        end
        function dimn = nz(obj)
            %nzp: number of nonlinearoutputs
            dimn = obj.P.nz;
        end
        function dimn = nw(obj)
            %nzp: number of performance outputs
            dimn = obj.P.nw;
        end
        function dimn = nwp(obj)
            %nzp: number of performance outputs
            dimn = obj.P.nwp;
        end
        function dimn = ny(obj)
            %nzp: number of performance outputs
            dimn = obj.P.nw;
        end
        function dimn = nu(obj)
            %nzp: number of performance outputs
            dimn = obj.P.nu;
        end

        function dimn = nxi(obj)
            %nxi: number of states in controller
            dimn = length(obj.K{1}.A);
        end

        %% get component functions
        function Scurr = get_P(obj, param)
            %GET_P get the plant at the current parameters
            % Pcurr = obj.weight_sum(obj.P, param.par);            

            Ak = zeros(obj.nxn, obj.nxn);
            Bk = zeros(obj.nxn, obj.nw + obj.nwp + obj.nu);
            Ck = zeros(obj.nz + obj.nzp + obj.ny, obj.nxn);
            Dk = zeros(obj.nz + obj.nzp + obj.ny, obj.nw + obj.nwp + obj.nu);
            for i = 1:obj.Nss   

                [Ac, Bc, Cc, Dc] = ssdata(obj.P{i});
                Ak = Ak + Ac * param.par(i);
                Bk = Bk + Bc * param.par(i);
                Ck = Ck + Cc * param.par(i);
                Dk = Dk + Dc * param.par(i);
                
            end
            Scurr = ss(Ak, Bk, Ck, Dk, 1);
            % Scurr = Pcurr.ss();
        end

        function [A, B1, B2, C1, D11, D12, C2, D21, D22] = ss_zy_wu(obj, param)
            %fetch matrices from each subsystem
            
            A = cell(obj.Nss, 1);
            B1 = cell(obj.Nss, 1);
            B2 = cell(obj.Nss, 1);
            C1 = cell(obj.Nss, 1);
            D11 = cell(obj.Nss, 1);
            D12 = cell(obj.Nss, 1);
            C2 = cell(obj.Nss, 1);
            D21 = cell(obj.Nss, 1);
            D22 = cell(obj.Nss, 1);
            for i = 1:obj.Nss
                [A{i}, B1{i}, B2{i}, C1{i}, D11{i}, D12{i}, C2{i}, D21{i}, D22{i}] = obj.P.ss_zy_wu(i);
            end

            A = obj.weight_sum(A, param);
            B1 = obj.weight_sum(B1, param);
            B2 = obj.weight_sum(B2, param);
            C1 = obj.weight_sum(C1, param);
            D11 = obj.weight_sum(D11, param);
            D12 = obj.weight_sum(D12, param);
            C2 = obj.weight_sum(C2, param);
            D21 = obj.weight_sum(D21, param);
            D22 = obj.weight_sum(D22, param);

        end

        function corn = interp(obj, par_desired)
            %vertex interpolator

            Npar = obj.Nss;
            Ncorn = obj.Ncorner;

            %solve using linear programming

            Aeq = [ones(1, Ncorn); obj.par'];
            beq = [1, par_desired];

            Aineq = -eye(Ncorn);
            bineq = zeros(Ncorn, 1);


            
            cost = zeros(Ncorn, 1);


            corn = linprog(cost, Aineq, bineq, Aeq, beq);
            % lb = zeros(Ncorn, 1);
            % ub = Inf*ones(Ncorn, 1);
            % Aineq = eye(Ncorn);
            % bineq = 

        end

        function Scurr = get_K(obj, param)
            %GET_K get the controller at the current parameters

            corn = obj.interp(param.par);

            
            % K_block = zeros(obj.nxi + obj.ny, obj.nxi + obj.nu);
            Ak = zeros(obj.nxi, obj.nxi);
            Bk = zeros(obj.nxi, obj.ny);
            Ck = zeros(obj.nu, obj.nxi);
            Dk = zeros(obj.nu, obj.ny);
            for i = 1:obj.Ncorner   

                [Ac, Bc, Cc, Dc] = ssdata(obj.K{i});
                Ak = Ak + Ac * corn(i);
                Bk = Bk + Bc * corn(i);
                Ck = Ck + Cc * corn(i);
                Dk = Dk + Dc * corn(i);
                % curr_block = corn(i) * [Ac, Bc; Cc, Dc];
            end
            
            Scurr = ss(Ak, Bk, Ck, Dk, 1);            
        end  

        %% accumulators
        function P_accum = weight_sum(obj, P, parl)
            %weighted sum over the parameters
            
            P_accum = zeros(size(P{1}));
            for i = 1:length(parl)
                P_accum = P_accum + parl(i) * P{i};
            end
        end


        function P_corn = weight_sum_corners(obj, P)
            %weighted sum over each corner of the polytope
            P_corn = cell(obj.Ncorner, 1);
            for i = 1:obj.Ncorner   
                P_corn = obj.weight_sum(P, obj.par(i, :));
            end
        end        

        %% tracking
        function [Sbeta, Rbeta] = get_tracked_opt(obj, param)
            %GET_TRACKED_OPT get the tracked position of the optimal
            %solution. allow for time-varying exosystems (periodic),
            %represented by a cell

            if isempty(obj.tracking)
                Sbeta0 = 1;
                Rbeta0 = 1;
            else
                Sbeta0 = obj.tracking.Sbeta;
                Rbeta0 = obj.tracking.Rbeta;
            end
 
            if nargin < 2 || isempty(param)
                if iscell(Sbeta0)
                    Sbeta = Sbeta0;
                    Rbeta = Rbeta0;
                else
                    Sbeta = cell(obj.Ncorner, 1);
                    Rbeta = cell(obj.Ncorner, 1);
    
                    
                    for i = 1:obj.Ncorner
                        Sbeta{i} = obj.weight_sum_corners(Sbeta0);
                        Rbeta{i} = obj.weight_sum_corners(Rbeta0);
                    end
                end
            else
                Sbeta = obj.weight_sum(Sbeta0, param);                
                Rbeta = obj.weight_sum(Rbeta0, param);                
            end
        end

        %% build the plant
        function [alg_psi, iqc_op, alg_loop] = build_plant(obj, iqc_data, rho)
            %BUILD_PLANT: form the plant to be used for analysis
            %or synthesis
            %Input:
            %   iqc_data: from manager.iqc_op_all, information about the
            %             operator iqc descriptions
            %   task: 'analysis' or 'synthesis'
            %
            %Output:
            %   alg_psi:    plant with filters (psi)
            %   alg_loop:   plant without filters, but after loop
            %               transformation (should be stable)
            %   iqc_op:     iqcs for the robust uncertainties


            alg_psi = cell(obj.Ncorner, 1);
            alg_loop = cell(obj.Ncorner, 1);

            %get the plant and the IQCs.
            %for each subsystem
            %
            %TODO: different IQCs for each subsystem, right now they are
            %identical.
            %

            ds = obj.get_discount();
            for i = 1:obj.Ncorner

                rho_curr = rho;

                param = struct('par', obj.par(i, :));
                alg = obj.get_alg(param); 
                [alg_psi{i}, iqc_op, alg_loop{i}] = build_plant_single(obj, alg, iqc_data, rho_curr);
            end
            %repeat this call multiple times for switched systems. This
            %function will be overloaded, whereas build_plant_single will
            %stay the same.




        end

        function tp = get_type(obj)
            %get the type of the switched system
            %is lpv_poly!
            tp = obj.type;
        end

        function mode_next = next_mode(obj, mode)
            %next mode in the switching sequence
            mode_next = [];
            error('need to implement the next mode nicely')
        end
    end
end


