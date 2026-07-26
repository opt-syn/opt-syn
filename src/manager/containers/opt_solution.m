classdef opt_solution
    %OPT_SOLUTION solution to an Analysis or Synthesis program
    %acquired from running an opt_manager object


    properties
        info; %information about the LMI solution
        dia; %feasibility of solution, should be negative
        objective; %minimization objective
        status; %zero if feasible, nonzero if infeasible
        lmi_out; %LMI constraints in the solution
        sys; %final algorithmic interconnection
        rho; %certified linear convergence rate
        cert=[]; %other certificates of the solution
        vars; %recovered variables of the problem
        recovery = []; %optional recovery of LMI constraints, blocks and their minimal eigenvalues
    end

    methods
        function obj = opt_solution(task)
            %OPT_SOLUTION constructor   
 
            if nargin 
                if  strcmp(task, 'analysis')
                    obj.cert = cert_analysis;
                else
                    obj.cert = cert_synthesis;
                end
            else
                obj.cert = struct;
            end
        end
    end
end