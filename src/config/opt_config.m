classdef opt_config
    %OPT_CONFIG configuration options for optimization algorithm analysis
    %and synthesis
    
    %these are separated into classes to allow for autodoc.
    properties       
        gen = opt_config_gen(); %general options for analysis and synthesis

        %analysis only options        
        ana = opt_config_ana(); %analysis options

        %synthesis only options
        syn = opt_config_syn(); %synthesis options

        
        % syn = struct('reduced_order', true, ... 
        %     'D_mask', [], ...           %which elements of D can be nonzero?
        %     'elimination', false,...);  %if there is only one specification 
        %     ...                         % (and the system is LTI), use the 
        %     ...                         % matrix elimination lemma    
        %     'elimination_type', 2);     %the elimination type                                        
        %                                 % 2: remove [Ak, Bk; Ck, Dk]
        %                                 % 1: remove [Ak, Bk; Ck1, Dk1]
        %                                 % 0: remove [Ak; Ck]
        %                                 %
        %                                 % true: smaller size variables
        %                                 % false: more interpretable

        %numerical tolerances for solving LMIs
        tol = opt_config_tol();             %upper bound on norm of controller state space matrices


        %recovery of the solution
        recovery = struct('blocks', false) %recover LMI blocks of the solution?

        %system-type specific options
        switched = struct('common', false); %use a common lyapunov function?

        %bisection options
        bisect = opt_config_bisect(); 
    end
    
    methods
        function obj = opt_config()
            %OPT_CONFIG construct the configuration            
        end
        
        function verdict = LMILAB(obj)
            %is LMILAB used?
            verdict = strcmp(obj.gen.solver, 'lmilab');
        end
    end
end

