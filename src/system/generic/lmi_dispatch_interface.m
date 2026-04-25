classdef lmi_dispatch_interface
    %LMI_DISPATCH_INTERFACE analysis and synthesis LMIs for the
    %algorithmic interconnections
    %
    %This contains generic routines common among both analysis and
    %synthesis for every system type
    %  
    properties
        sys;
        LMILAB = 1;
        tol = struct('M', 1e-7, ... %tolerance for dissipation constraints
            'X', 1e-7, ...          %tolerance for sign/terminal cost constraints 
            'G_max', 100)           %upper bound on norm of storage matrix
    end
    
    methods
        function obj = lmi_dispatch_interface(sys)
            %LMI_DISPATCH_INTERFACE Construct an instance of this class
            %   Detailed explanation goes here
            obj.sys = sys;
        end
        
        
        %KYP lemma terms, commonly found matrices
        
        function [cons, objective, con_M] = con_dissipation(obj, vars, cons, diss)
            %CON_DISSIPATION form the dissipation constraint
            %
            %Input:
            %   vars:   variables of the problem        
            %   cons:   accumulated constraints
            %   diss:   structure describing the problem
            %       plant:  system to control
            %       spec:   performance specification           
            %       target: whether the performance measure should be optimized
            %               true:  soft constraint (e.g. Schur complement
            %                                       formulation)
            %               false: hard constraint            
            %
            %Output:
            %   cons:   accumulated constraints
            %   objective:  term to be minimized            
            
            
            %need to look up the right constraint            

            %Upper-levels: iterate over the systems
            [cons, objective, con_M] = obj.con_dissipation_single(vars, cons, diss);

            
                      
        end
    

        function [cons, objective, con_M] = con_dissipation_single(obj,  vars, cons, diss)
            %CON_DISSIPATION form a single dissipation constraint
            %
            %Input:
            %   vars:   variables of the problem            
            %   diss:   structure describing the problem
            %       plant:  system to control
            %       spec:   performance specification           
            %       target: whether the performance measure should be optimized
            %               true:  soft constraint (e.g. Schur complement
            %                                       formulation)
            %               false: hard constraint
            %   param:  other parameters
            %
            %Output:        
            %   objective:  term to be minimized            
    
    
            %need to look up the right constraint            
    
            lmi_handle = str2func(diss.spec.type);
    
            [cons, objective, con_M] = lmi_handle(obj, vars, cons, diss);
    
        end

        %% helper functions
        function M = merge_spec_M(obj, iqc_rob, sp)
            %MERGE_SPEC_M merge the running cost of the robustness and the 
            %performance specification
            %
            %Input:
            %   iqc_rob: robust IQC 
            %   sp:      performance specification

            n1 = iqc_rob.np;
            m1 = iqc_rob.nq;
            n2 = length(sp.iwp);
            m2 = length(sp.izp);
            [M] = outer_blkdiag(iqc_rob.M, sp.supply, n1, m1, n2, m2);
            % Mdiag = blkdiag(iqc_op.M, sp.supply);
        end
        
    end




    methods (Abstract)
        create_vars(obj, vars, cons, alg, specs)
        
        
        %supported operators
        stability(obj, vars, cons, diss)
        
        
    end
end

