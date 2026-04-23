classdef opt_spec_inteface
    %OPT_SPEC performance specification for a system
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
        vars;
        type;
        bound;
        LMILAB = 1;        
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

            if strcmp(type, 'stab') || strcmp(type, 'stability') || strcmp(type, 'finite_l2')
                obj.rho = bound;
            end

        end

        function [vars, cons] = create_vars(obj, cons)
        end
        
        function [vars, cons, iqc] = create_iqc(obj, cons)
            %CREATE_IQC Summary of this method goes here
            %   Detailed explanation goes here

            if nargin < 2;
                cons = [];
            end
            nwp = length(obj.iwp);
            nzp = length(obj.izp);
            switch obj.type
                case 'e2e'
                    vars = []; iqc = iqc_e2e(nwp, nzp, obj.bound);
                case 'finite_l2'
                    mu_l2 = lmim('mu_l2', 1, 1);
                    
                    cons = append_lmi(cons, mu_l2, obj.LMILAB);

                    % cons = append_lmi(cons, mu_l2, obj.LMILAB);
                    cons = append_lmi(cons, obj.finite_l2_bound - mu_l2, obj.LMILAB);

                    vars = struct('mu_l2', mu_l2); 
                    % iqc = iqc_finite_l2(nwp, nzp, mu_l2);
                    % iqc = iqc_finite_l2(nwp, nzp, 2500);
                    % iqc = iqc_e2e(nwp, nzp, 50);
                    iqc = iqc_finite_l2(nwp, nzp, 50);
                otherwise
                    vars = []; iqc = [];
            end
        end
    end
end

