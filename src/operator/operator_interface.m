classdef (Abstract) operator_interface
    %OPERATOR_INTERFACE An operator with relation w \in F(z)
    % container for the variables, IQCs, and relations    
    
    properties
        % vars; %variables in the analysis/synthesis programs
        id=1;   %the ID of the operator in the problem (index)
        % iqc;  %IQC descring relations satisfied by the operator
        %monotone = true; %can the operator be bounded using monotonicity?

        c = 1; %coordinate dimension for the operator
        EQUALITY = false;
        LMILAB = true;
    end
    
    methods
        function obj = operator_interface(c)
            %OPERATOR_INTERFACE Construct an instance of this class
            %   Detailed explanation goes here
            if nargin  > 0
                
                obj.c = c;
            end
        end 

        function strid = sid(obj)
            strid = num2str(obj.id);
        end

        function [iqc, vars, cons] = create_iqc(obj, cons, order, reps)
            %CREATE_VARS form the IQC for the general operator

            %Input: 
            %   order:  order of the IQC [number of lags]
            %   rep:    number of repetitions of the operator (non-frugal)
            %
            %Output:
            %   vars:   variables of the problem
            %   cons:   constraints in the problem (in terms of the
            %           variables directly)

            if obj.same
                iqc = obj.get_same(reps);
                vars = [];
            else
                [vars] = obj.create_vars(order, reps);
      
                %form the IQC            
                M = build_M(obj, vars, order, reps);
                X = build_X(obj, vars, order, reps);
    
                loop = obj.build_loop(reps);            
    
                [psi1, psi2] = obj.build_psi(vars, order, reps);
    
                iqc_orig = iqc_loop_split(psi1, M, loop, psi2, X);
    
                iqc = iqc_orig.lift(obj.c);
                cons = obj.filter_constraints(cons, order, vars, iqc);
            end
        end

        function sm = same(obj)
            %SAME: is there any uncertainty in this oracle?
            %Overriden by class Sml
            sm = false;
        end

        function sm = get_same(obj)
            %GET_SAME: is there any uncertainty in this oracle?
            %Overriden by class Sml
            sm = [];
        end
    end

    methods (Abstract)
        %constructor
        create_vars(obj, order, reps);        
        
        build_psi(obj, vars, order, reps);
        build_M(obj, vars, reps);
        build_X(obj, vars, order, reps);
        build_loop(obj, reps);

        filter_constraints(obj, cons, order, vars, iqc)
        %recovery
%         factor_iqc(obj)        
    end

end