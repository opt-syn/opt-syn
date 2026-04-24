%OLD CODE OLD CODE OLD CODE OLD CODE
%
%This should be used as a reference for future implementation
%
classdef alg_options
    %ALG_PROBLEM configuration options for algorithm analysis and synthesis
    
    %
    %
    properties
        %problem definitions
        m = [];
        L = [];
        gen = false; %used modified code (for general iqc)

        %common options
        solver = 'LMILAB';
        verbose = 0;
        tol = 1e-4;
        pass_tol = 0.001;

        L_top = []; %should f_L be filtered (default, = 1), or f_m be filtered ( = 0)

        rho = 1;
        rho_range = [1e-4, 1.01];
        % analysis
        mul_order = 0;
        mul_type = 1;

        repeated = []; %cell of indices for repeated nonlinearities
        full_block = false;

        common = true;

        % synthesis
        static = 0;
        D_mask = []; %sparsity pattern on D controller
        recover = true;
        HINF = false;
        spread_tol = 0.01;        
        elim = 0;
        gam_max = 500;
        mul = {};
        minreal_tol = 1e-5;
        onesign = -1; %sign in front of the constant offset
        
        %static subcontroller initial point
        b_start = [];

        %internal model configuration
        tracking = [];
        basis_trans = false; %coordinate transformation

    end
    
   
end

