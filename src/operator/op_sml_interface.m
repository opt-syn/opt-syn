classdef op_sml_interface < operator_interface
    %OP_SML An operator which is the subdifferential of a function in SmL:
    
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
        m; %lower bound parameter
        L; %upper bound parameter
        L_top = true; %should the L term be filtered (Lz - w)?
        ERGODIC = false; %function value penalties?
    end
    
    methods
        function obj = op_sml_interface(m, L, c)
            %OP_SML_INTERFACE Constructor
            %
            %Args:
            %   m: lower bound parameter
            %   L: upper bound parameter
            %   c: dimension of coordinate lift

            if nargin < 3
                c = 1;
            end
            obj@operator_interface(c)            

            obj.m = m;
            obj.L = L;
        end
       

        % common routines

        function se = same(obj)
            %SAME is there any uncertainty in this oracle?            
            %
            %Returns:
            %   sm (bool): m=L?.
            se = obj.L == obj.m;
        end
        
        function sig = sigma(obj)            
            %SIGMA used to define all IQCs
            %
            %Return:
            %   sig: 1/(L-m)
            sig = 1/(obj.L - obj.m);
        end
               

        function loop = build_loop(obj, reps)
            %BUILD_LOOP construct the signal transformation matrix
            %
            %Args:
            %   reps:    number of repetitions of the operator (from the bind)
            %
            %Returns:
            %   loop_out: signal transformation matrix for the operator

            if obj.L_top
                loop_base = [-obj.sigma, 1; 1, obj.m];
            else
                loop_base = [-obj.sigma, 1; -1, obj.L];
            end

            loop = kron(loop_base, eye(reps));
        end


        function M_erg = ergodic_supply(obj, reps)
            %ERGODIC_SUPPLY supply rate for function value decrease
            %ergodic convergence             
            %
            %Args:
            %   reps:    number of repetitions of the operator (from the bind)
            %
            %Returns:
            %   M_erg: quadratic supply for function values

            
            if nargin < 2
                reps = 1;
            end
            if obj.same
                M_erg = [];
            else
                sig = 1/(obj.L - obj.m);

                m = obj.m;

                mp = m*[1, sig; sig, sig^2];
                mq = sig*[0, 0; 0, 1];
                M_erg = kron(mp + mq, eye(reps));
            end
        end


        function [iqc] = create_iqc_identity(obj, reps)
            %CREATE_IQC_IDENTITY form a valid IQC satisfied by the general
            %operator. This is used as a warm start in synthesis.
            %
            %Args:             
            %   reps:    number of repetitions of the operator (from the bind)
            %
            %Returns:
            %   iqc (iqc_loop_split): a valid IQC with no dynamics    

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



        function sm = get_same(obj, reps)
            %GET_SAME explicit matrix in LFT           
            %
            %Returns:
            %   sm: m==L? If so, return the explicit matrix mI.
            sm = kron(obj.m, eye(reps));
        end
        
    end

    
end

