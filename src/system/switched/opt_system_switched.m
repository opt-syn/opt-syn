classdef  opt_system_switched < opt_system_interface
    %OPT_SYSTEM_SWITCHED interconnection of network and operators
    %polytopic setting: a cell A = sum theta_i A_i for parameters theta_i
    %
    %useful for switched systems, periodic systems, and LPV systems
    
    
    properties        
        adj;    %switching graph  (adjacency matrix)              
    end
    
    methods
        function obj = opt_system_switched(op, P, K, adj, bind, tracking)
            %OPT_SYSTEM constructor            
            if nargin < 5
                s = length(op);
                bind = 1:s;            
            end

            if nargin < 6
                tracking = [];
            end

            if ~iscell(K)
                K0 = K;
                K = cell(P.Nss, 1);
                for i = 1:P.Nss
                    K{i} = K0;
                end
            end

            if nargin >= 5
                 tracking = [];
            end

            obj@opt_system_interface(op, P, K, bind, tracking)            
            obj.adj = adj;
            obj.type = 'switched';
        end        

        %TODO: allow for parameterized systems

        function [src, dst] = get_arcs(obj)
            %GET_ARCS get transitions in the adjacency matrix
            [src, dst] = find(obj.adj);
        end

        function tp = get_type(obj)
            %get the type of the switched system
            if all((obj.adj==0) + (obj.adj==1), 'all')
                %robust switching
                %TODO: not yet implemented
                tp = 'switched';
            else
                %TODO: not yet implemented
                %stochastic: markov jump linear system
                tp = 'mjls';
            end   
        end

        function mode_next = next_mode(obj, mode)
            if nargin < 2
                mode = 1;
            end

            %TODO: debug this
            g = obj.adj(mode, :);
            gc = g/sum(g);

            gs = cumsum(gc);

            u = rand(1);

            % mode_next = [];
            mode_next = find(u <= gs, 1, 'first');
        end

        function Scurr = get_P(obj, param)
            %GET_P get the plant
            if isnumeric(param) 
                if isempty(param)
                    Scurr = obj.P;
                else
                    Scurr = obj.P.P{param};
                end
            else
                Pcurr = obj.P{param.mode};
                Scurr = Pcurr.ss();
            end
        end

        function Kcurr = get_K(obj, param)
            %TODO: override this with parameters
            if isnumeric(param ) 
                if isempty(param)                    
                    Kcurr = obj.K;
                else
                    Kcurr = obj.K{param};
                end
            else
                if isa(obj.K{param.mode}, 'genplant')
                    Kcurr = obj.K{param.mode}.ss;
                else
                    Kcurr = obj.K{param.mode};
                end
            end
        end    
        
        %dimensions
        function dimn = nxn(obj)
            %nxn: number of states in network
            dimn = obj.P.nx;
        end

        function ds = get_discount(obj)
            %which subsystems are exponentially discounted?
            ds = obj.discount;
            if length(ds)==1                
                ds = ones(1, obj.Nss)*ds;
            end
        end

        function nss = Nss(obj)
            nss = obj.P.Nss;
        end

        function dimn = nxi(obj)
            %nxi: number of states in controller
            dimn = length(obj.K{1}.A);
        end

        function pow = discount_schedule(obj, ordermax)
            %DISCOUNT_SCHEDULE exponential weights encountered when
            %applying the FIR filters
            %Args:
            %   ordermax: maximum order of the IQCs
            %
            %Return:
            %   pow: Exponent sequence of discounts
            %
            % Example: 
            %   [0; 1; 2] -> rho.^[0; 1; 2] for uniform exponential stability
            %   [0; 0; 1] -> rho.^[0; 0; 1] for shuffled switched stability
            %
            
            %This becomes relevant when performing shuffled systems
            %(override on switched systems) 
            


            if isscalar(obj.discount)
                pow = -(0:ordermax)';
            else


                Gcon = obj.adj;
                u_discount = [];
                for i = 1:obj.Nss

                    all_walks_curr= all_walks(Gcon, ordermax, i, []);


                    discount_curr = obj.discount(all_walks_curr);
                    u_curr = unique(discount_curr, "rows");

                    u_discount = unique([u_discount; u_curr], "rows");

                    % all_walks = [all_walks, all_walks_curr]
                end

                pow = -cumsum(u_discount, 2)';
            end


        end


        %% build the plant
        function [alg_psi, iqc_op, alg_loop] = build_plant(obj, iqc_data)
            %BUILD_PLANT: form the plant to be used for analysis or synthesis
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


            alg_psi = cell(obj.Nss, 1);
            alg_loop = cell(obj.Nss, 1);
           
            %get the plant and the IQCs.
            %for each subsystem
            %
            %TODO: different IQCs for each subsystem, right now they are
            %identical.
            %

            ds = obj.get_discount();
            for i = 1:obj.Nss

 

                param = struct('mode', i);
                alg = obj.get_alg(param); 
                [alg_psi{i}, iqc_op, alg_loop{i}] = build_plant_single(obj, alg, iqc_data);
            end
            %repeat this call multiple times for switched systems. This
            %function will be overloaded, whereas build_plant_single will
            %stay the same.
            

                      

        end              

    end
end

