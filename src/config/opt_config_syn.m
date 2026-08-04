classdef opt_config_syn
    %OPT_CONFIG_SYN Configuration options for synthesis
    
    properties
        reduced_order = true; %use internal model structure to synthesize reduced-order controllers (LTI and periodic-orbit)               
        D_mask = []; %which elements of D can be nonzero?
        elimination = true; % if there is only one specification 
             % (and the system is LTI), use the 
            % matrix elimination lemma
        elimination_type = 2; %the elimination type (no longer required)
    end
    
    methods
        function obj = opt_config_syn()
            %OPT_CONFIG_SYN constructor
        end        
    end
end

