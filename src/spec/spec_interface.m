classdef spec_interface
    %SPEC performance specification for an optimization algorithm system    
    %
    %default: stability
    %others include:
    %   energy to energy gain (ell 2)
    %   energy to peak gain   (generalized H 2)
    %   peak to peak gain   (generalized H 2)
    %   Covariance amplification (H 2)

    %none of these specifications will involve loop transformations
    %   (simplification of the routines)
    
    properties
        rho=1; %exponential discounting
        iwp = [];    %indices for performance input
        izp = [];    %indices for performance output        
        vars;                
        LMILAB = 1;            
        id = 0;
        target = false; %should this be the target of optimization (within 
                        %bisection iterations)?
    end
    
    methods
        function obj = spec_interface(iwp, izp, rho)
            %OPT_PERFORMANCE a performance specification for the IQC
            %analysis.
            %
            %Input:
            %   type:   the specification (e.g. 'stability', 'e2e')
            %   bound:  the current value of the bound
            %   iwp:    performance inputs in the network    
            %   izp:    performance outputs in the network            

            if nargin > 1
                obj.iwp = iwp;
                obj.izp = izp;
            end
            if nargin > 2
                obj.rho = rho;
            end
            
            
        end


        function nzzp = nzp(obj)
            %number of outputs
            nzzp = length(obj.izp);
        end

        function nwwp = nwp(obj)
            %number of inputs
            nwwp = length(obj.iwp);
        end

        function [vars, cons] = create_vars(obj, cons, name)
            %CREATE_VARS form the variables for the problem            
            vars = [];            
        end

        function [M] = supply(obj, vars_spec)
            %SUPPLY quadratic performance specification
            M = [];
        end

        function [obj] = set_p(obj, p)
            %SET_P set a parameter when performing bisection
            %
            %
            %Example: Peak-to-Peak norm certifier
            %or l2 Gain bound
            
        end

    end

    methods (Abstract)
        % create_vars(obj, cons)
        % supply(obj)
    end
        
        % function [vars, cons, iqc] = create_iqc(obj, cons)
        %     %CREATE_IQC Summary of this method goes here
        %     %   Detailed explanation goes here
        % 
        %     if nargin < 2;
        %         cons = [];
        %     end
        %     nwp = length(obj.iwp);
        %     nzp = length(obj.izp);
        %     switch obj.type
        %         case 'e2e'
        %             vars = []; iqc = iqc_e2e(nwp, nzp, obj.bound);
        %         case 'finite_l2'
        %             mu_l2 = lmim('mu_l2', 1, 1);
        % 
        %             cons = append_lmi(cons, mu_l2, obj.LMILAB);
        % 
        %             % cons = append_lmi(cons, mu_l2, obj.LMILAB);
        %             cons = append_lmi(cons, obj.finite_l2_bound - mu_l2, obj.LMILAB);
        % 
        %             vars = struct('mu_l2', mu_l2); 
        %             % iqc = iqc_finite_l2(nwp, nzp, mu_l2);
        %             % iqc = iqc_finite_l2(nwp, nzp, 2500);
        %             % iqc = iqc_e2e(nwp, nzp, 50);
        %             iqc = iqc_finite_l2(nwp, nzp, 50);
        %         otherwise
        %             vars = []; iqc = [];
        %     end
        % end
    
end

