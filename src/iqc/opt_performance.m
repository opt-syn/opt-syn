classdef opt_performance
    %OPT_PERFORMANCE performance specification for a system
    %
    %default: stability
    %others include:
    %   energy to energy gain (ell 2)
    %   energy to peak gain   (generalized H 2)
    %   peak to peak gain   (generalized H 2)
    %   Covariance amplification (H 2)
    
    properties
        rho=1; %exponential discounting
        iwp = [];    %indices for performance input
        izp = [];    %indices for performance output
        iqc;
        vars;
        type;
        bound;
    end
    
    methods
        function obj = opt_performance(type, bound, iwp, izp)
            %OPT_PERFORMANCE a performance specification for the IQC
            %analysis.
            %
            %Input:
            %   type:   the specification (e.g. 'stability', 'e2e')
            %   bound:  the current value of the bound
            %   iwp:    performance inputs in the network    
            %   izp:    performance outputs in the network


            obj.type = type;

            if nargin < 2
                bound = 1;
            end
            if nargin > 3
                obj.iwp = iwp;
                obj.izp = izp;
            end
            obj.bound = bound;

            if strcmp(type, 'stab') || strcmp(type, 'stability')
                obj.rho = bound;
            else
                [obj.vars, obj.iqc] = obj.create_iqc();
            end
            
            


        end
        
        function [vars, iqc] = create_iqc(obj)
            %CREATE_IQC Summary of this method goes here
            %   Detailed explanation goes here

            nwp = length(obj.iwp);
            nzp = length(obj.izp);
            switch obj.type
                case 'e2e'
                    vars = []; iqc = iqc_e2e(nwp, nzp, obj.bound^2);

            end
        end
    end
end

