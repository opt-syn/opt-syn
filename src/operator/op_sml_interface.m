classdef op_sml_interface < operator_interface
    %OP_SML An operator which is the subdifferential of a function in SmL:
    %
    %F = partial f, where f(x) - m norm(x, 2)^2 and L norm(x, 2)^2 - f(x)
    %are both proper, convex, and closed with -Inf < m <= L <= Inf
    %
    % if m=0, L=inf: convex
    % if m>0:        m-strongly convex
    % if m<0:        (-m)-weakly convex
    % if m>0, L<inf: m-strongly convex and L-smooth (L-Lipschitz gradients)
    %
    %
    % TODO: generalize to matrices m and L?
    
    properties
        m; %convexity parameter
        L; %smoothness parameter
        L_top = true; %should the L term be filtered (Lz - w)?
    end
    
    methods
        function obj = op_sml_interface(m, L, c)
            %OP_SML Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 3
                c = 1;
            end
            obj@operator_interface(c)            

            obj.m = m;
            obj.L = L;
        end
       

        % common routines

        function se = same(obj)
            %SAME no loop transformation or IQC required
            %perfectly known oracle
            se = obj.L == obj.m;
        end
        
        function sig = sigma(obj)            
            %SIGMA used to define all IQCs
            sig = 1/(obj.L - obj.m);
        end
               

        function loop = build_loop(obj, reps)
            %BUILD_LOOP construct the loop transformation
            if obj.L_top
                loop_base = [-obj.sigma, 1; 1, obj.m];
            else
                loop_base = [-obj.sigma, 1; -1, obj.L];
            end

            loop = kron(loop_base, eye(reps));
        end


        function [iqc] = create_iqc_identity(obj, reps)
            %CREATE_VARS form the IQC for the general operator
            %identity IQC in psi

            %Input:             
            %   rep:    number of repetitions of the operator (non-frugal)
            %
            %Output:
            %   vars:   variables of the problem
            %   cons:   constraints in the problem (in terms of the
            %           variables directly)

            if nargin < 2
                reps = 1;
            end
            
            if obj.same
                iqc = obj.get_same(reps);
            else
                vars = [];
                %form the IQC            
                % M = build_M(obj, vars, 0, reps);
                M = kron([0, 1; 1, 0], eye(reps));
                X = [];
    
                loop = obj.build_loop(reps);            
    
                % [psi1, psi2] = obj.build_psi(vars, order, reps);
                psi1 = ss(eye(reps));
                psi2 = ss(eye(reps));
    
                iqc_orig = iqc_loop_split(psi1, M, loop, psi2, X);
    
                iqc = iqc_orig.lift(obj.c);                
            end
        end
        
    end
end

