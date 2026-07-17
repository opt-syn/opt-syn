classdef op_sml_causal_orbit < op_sml_causal
%OP_SML_CAUSAL_ORBIT 

%OP_SML_CAUSAL_ORBIT An operator which is the subdifferential of a function in SmL:
%
%F = partial f, where f(x) - m norm(x, 2)^2 and L norm(x, 2)^2 - f(x)
%are both proper, convex, and closed with -Inf < m <= L <= Inf
%
%
%Time-variation is included in the relation: 
%   w'_k \in R^(-k) \partial f(R^k z'_k)
%
%R is an orthogonal matrix

    properties 
        R = 1; %coordinate shift matrix
    end

    methods
        function obj = op_sml_causal_orbit(m, L, R)
            %OP_SML Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 3
                R = 1;
            end
            c = size(R, 1);
            obj@op_sml_causal(m ,L, c)            
    
            obj.R = R;
        end
        


        function [iqc, vars, cons] = create_iqc(obj, cons, order, reps)
            %CREATE_VARS form the IQC for the general operator

            %Input: 
            %   order:  order of the IQC [number of lags]
            %   reps:    number of repetitions of the operator (non-frugal)            
            %
            %Output:
            %   vars:   variables of the problem
            %   cons:   constraints in the problem (in terms of the
            %           variables directly)

            [iqc, vars, cons] = create_iqc@op_sml_causal(obj, cons, order, reps);

            n = size(iqc.Psi1.A, 1);
            Rkroninv = kron(eye(n/obj.c), obj.R');
            iqc.Psi1.A = Rkroninv * iqc.Psi1.A;
            iqc.Psi1.B = Rkroninv * iqc.Psi1.B;


        end
    end
end