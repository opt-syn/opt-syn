classdef regulator_interface
    %REGULATOR_INTERFACE 
    %
    % A regulator employed for the synthesis of optimization algorithms
    %
    % This regulator carries the necessary internal model required for
    % convergence of optimization algorithms (constant shift of optimal 
    % solution), as well as extra states associated with position of the
    % unknown disturbance. This regulator is adjoined to the network
    % dynamics in synthesis of algorithms.
    %
    % The regulator allows for the perfect tracking of optimal solutions
    % (if known) by the tracking dynamics (Sbeta, Rbeta) in (sys).
    %
    %
    % The regulator is specialized for a specific kind of system
    %
    % 
    %
    
    properties
        sys;

        %the exosystem for the optimal solution and operators
        S;
        R;

        %nominal solution to the regulator equations
        Pi;
        Gam;
        Phi;

        %freedom in solving the regulator equations
        Pi_basis;
        Gam_basis;
        Phi_basis;
    end
    
    methods
        function obj = regulator_interface(sys)
            %REGULATOR_INTERFACE Construct an instance of this class
            %   Detailed explanation goes here
            obj.sys = sys;

            obj = obj.form_internal_model();
        end

        function NS = ns(obj)
            %NS: number of states of exosystem

            NS = length(obj.S);
        end

        function vars_reg = create_vars(obj)
            vars_reg.Pi = obj.Pi;
            vars_reg.Gam = obj.Gam;
            vars_reg.Phi = obj.Phi;

            

            nbasis = size(obj.Gam_basis, 2);
            % if nbasis> 0 
            if false


                %search over solutions
                eta = lmim('reg_param', nbasis, 1, 'full');

                %TODO: parameterization of regulator equation solutions
                %needs to be debugged (e.g. Laplacian example)
                Pi_free_vec = obj.Pi_basis*eta;
                Gam_free_vec = obj.Gam_basis*eta;
                Phi_free_vec = obj.Phi_basis*eta;

                Pi_free = [];
                Gam_free = [];
                Phi_free = [];

                d = size(vars_reg.Pi, 2);
                n = size(vars_reg.Pi, 1);
                nu = size(vars_reg.Gam, 1);
                ny = size(vars_reg.Phi, 1);

                %parameterization of the nullspace
                for i = 1:d
                    Pi_curr = lmim_index(Pi_free_vec, (i-1)*n + (1:n), 1);
                    Gam_curr = lmim_index(Gam_free_vec, (i-1)*nu + (1:nu), 1);
                    Phi_curr = lmim_index(Phi_free_vec, (i-1)*ny + (1:ny), 1);

                    Pi_free = [Pi_free, Pi_curr];
                    Gam_free = [Gam_free, Gam_curr];
                    Phi_free = [Phi_free, Phi_curr];
                end

                %add the free nullspace part to the system

                vars_reg.Pi = vars_reg.Pi + Pi_free;
                vars_reg.Gam = vars_reg.Gam + Gam_free;
                vars_reg.Phi = vars_reg.Phi + Phi_free;
                vars_reg.reg_param = eta;
            else
                vars_reg.reg_param= [];
            end
        end
        
    end

    methods(Abstract)
        form_internal_model(obj)
        check_regulator(obj)
        get_model(obj)
    end
end

