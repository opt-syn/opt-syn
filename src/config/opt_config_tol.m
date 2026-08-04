classdef opt_config_tol
    %OPT_CONFIG_TOL Configuration options for numerical tolerances in analysis and synthesis
    
    properties
        dia= 1e-11;   %tolerance for acceptable solution 
        M= 1e-7;      %tolerance for dissipation constraints
        X=  1e-7;     %tolerance for sign/terminal cost constraints 
        G_max=  100;  %upper bound on norm of storage matrix (analysis)
        GX_max=  100; %upper bound on norm of primal storage matrix (synthesis)
        GY_max=  100; %upper bound on norm of dual storage matrix (synthesis)
        G_min=  1e-4; %lower-bound on storage function eigenvalue 
        spread=  0.01;%tolerance in the matrix dilation constraint
    input_diss=  1e-4;%tolerance for strict input dissipation
        K_max=  100; %upper bound on norm of controller state space matrices
    end
    
    methods
        function obj = opt_config_tol()
            %OPT_CONFIG_tol constructor
        end        
    end
end

