classdef  opt_system_lpv_poly < opt_system_switched
    %OPT_SYSTEM_LPV_POLY interconnection of network and operators
    %
    %an lpv-polytopic system: gain-scheduled controller

    
    properties
        %polytopes: convex hulls
        par = [];  %polytope of parameter values 
        pair = []; %polytope of transitions [parameter, next_paramater] 

        interp = []; %a vertex interpolator: returns a set of weights


        pair_h;

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

            %process the polytope set

            Npar = P.Nss;
            Ncorner = size(poly, 1);
            if size(ppar, 2) == Npar
                %all to all transitions
                %duplicate the polytope values
                obj.par = pair;
                obj.pair = [kron(ones(par, 1), poly), kron(par, ones(Npar, 1))];

                Npair = Npar * obj.par;
            else
                obj.ppair = ppair;
                par_iso = obj.ppair(:, 1:Npar);
                obj.par = unique(par_iso', 'rows')';
            end


            %vertex scheduled controller (default to the same controller
            %everywhere)
            if ~iscell(K)
                K0 = K;
                K = cell(Ncorner, 1);
                for i = 1:Ncorner
                    K{i} = K0;
                end
            end

            obj.interp = get_vertex_interp(obj.par);
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

        function dimn = nxi(obj)
            %nxi: number of states in controller
            dimn = length(obj.K{1}.A);
        end

        function Scurr = get_P(obj, param)
            %GET_P get the plant at the current parameters
            Pcurr = obj.weight_sum(P, param.par);            
            Scurr = Pcurr.ss();
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


        function P_corn = weight_sum_corners(obj)
            %weighted sum over each corner of the polytope
            P_corn = cell(obj.Ncorner, 1);
            for i = 1:obj.Ncorner   
                P_corn = obj.weight_sum(P, obj.par(i, :));
            end
        end        

        %% tracking
        function [Sbeta, Rbeta] = get_tracked_opt(obj)
            %GET_TRACKED_OPT get the tracked position of the optimal
            %solution. allow for time-varying exosystems (periodic),
            %represented by a cell

            [Sbeta0, Rbeta0] = get_tracked_opt@opt_system_switched(obj);

            if iscell(Sbeta0)
                Sbeta = Sbeta0;
                Rbeta = Rbeta0;
            else
                Sbeta = cell(obj.Ncorner, 1);
                Rbeta = cell(obj.Ncorner, 1);

                
                for i = 1:Ncorner
                    Sbeta{i} = obj.weight_sum_corners(Sbeta0);
                    Rbeta{i} = obj.weight_sum_corners(Rbeta0);
                end
            end
 
        end

        function tp = get_type(obj)
            %get the type of the switched system
            %is lpv_poly!
            tp = obj.type;
        end

        function mode_next = next_mode(obj, mode)
            %next mode in the switching sequence
            error('need to implement the next mode nicely')
        end
    end
end

