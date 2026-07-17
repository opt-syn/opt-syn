classdef  opt_system < opt_system_interface
    %OPT_SYSTEM interconnection of network and operators
    %by default is LTI (linear time invariant)
 
    
    methods
        function obj = opt_system(op, P, K, bind, tracking)
            %OPT_SYSTEM constructor for the system
            if nargin < 4
                bind = 1:length(op);
            end
            
            if nargin < 5
                 tracking = [];
            end

            obj@opt_system_interface(op, P, K, bind, tracking)
            
            obj.type = 'lti';
        end    

        %% Dimension Counters
        function dimn = nxn(obj)
            %nxn: number of states in network
            dimn = length(obj.P.A);
        end


        function dimn = nxi(obj)
            %nxi: number of states in controller
            dimn = length(obj.K.A);
        end


        

        %% getters
        function Pcurr = get_P(obj, param)
            %GET_P get the network P
            Pcurr = obj.P.ss();
        end
        
        function Kcurr = get_K(obj, param)
            %GET_K get the controller K
            if isa(obj.K, 'genplant')
                Kcurr = obj.K.ss;
            else
                Kcurr = obj.K;
            end
        end

        

        %% build the plant
        function [alg_psi, iqc_op, alg_loop] = build_plant(obj, iqc_data, rho)
            %BUILD_PLANT: form the plant to be used for analysis
            %or synthesis
            %Input:
            %   iqc_data: from manager.iqc_op_all, information about the
            %             operator iqc descriptions
            %   rho: exponential convergence rate (default 1)
            %
            %Output:
            %   alg_psi:    plant with filters (psi)
            %   alg_loop:   plant without filters, but after loop
            %               transformation (should be stable)
            %   iqc_op:     iqcs for the robust uncertainties

            if nargin < 3
                rho = 1;
            end

            

            %get the plant and the IQCs.
            
            if strcmp(iqc_data.task, 'analysis');
                alg = obj.get_alg([]);            
            else
                alg = obj.P.P;
            end
            

            %repeat this call multiple times for switched systems. This
            %function will be overloaded, whereas build_plant_single will
            %stay the same.
            [alg_psi, iqc_op, alg_loop] = build_plant_single(obj, alg, iqc_data, rho);

                      

        end 

        


    end
end

