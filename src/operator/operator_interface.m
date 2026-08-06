classdef (Abstract) operator_interface
    %OPERATOR_INTERFACE An operator with relation :math:`w \in F(z)`. 
    % A container for the variables, IQCs, and relations    
    
    properties        
        id=1;   %the ID of the operator in the problem (index)
        c = 1; %coordinate dimension for the operator
        EQUALITY = false; %is this an equality constraint?
        LMILAB = true; %is LMIlab used?
    end
    
    methods
        function obj = operator_interface(c)
            %OPERATOR_INTERFACE Constructor
            %
            %Args
            %   c (int): coordinate dimension for the operator
            if nargin  > 0
                
                obj.c = c;
            end
        end 

        function strid = sid(obj)
            %SID the convert the id to a string
            %
            %Return:
            %   strid (char): string of the id
            strid = num2str(obj.id);
        end

        function [iqc, vars, cons] = create_iqc(obj, cons, order, reps)
            %CREATE_IQC form the IQC for the operator class
            %
            %Args: 
            %   cons: accumulated constraints
            %   order:  order of the IQC [number of lags]
            %   reps:    number of repetitions of the operator (in the bind)
            %
            %Returns:
            %   iqc:    a valid iqc for the operator
            %   vars:   variables of the problem
            %   cons:   constraints in the problem (in terms of the
            %           variables directly)

            if obj.same
                iqc_orig = obj.get_same(reps);   
                iqc = kron(iqc_orig, eye(obj.c));
                vars = [];
            else
                if isscalar(order)
                    order = [order, 0];
                end
                

                [vars] = obj.create_vars(order, reps);
      
                %form the IQC            
                M = build_M(obj, vars, order, reps);
                X = build_X(obj, vars, order, reps);
    
                loop = obj.build_loop(reps);            
    
                [psi1, psi2] = obj.build_psi(vars, order, reps);
    
                iqc_orig = iqc_loop_split(psi1, M, loop, psi2, X);
    
                iqc = iqc_orig.lift(obj.c);
                
                % cons = obj.filter_constraints(cons, order, vars, iqc);
            end
            
        end

        function sm = same(obj)
            %SAME is there any uncertainty in this oracle?
            %Overriden by class Sml
            %
            %Returns:
            %   sm (bool): m=L? default to false.
            sm = false;
        end

        function sm = get_same(obj)
            %GET_SAME what is the explicitly known loop transformation matrix
            %Overriden by class Sml
            %
            %Returns:
            %   sm: m=L? default to empty.
            sm = [];
        end
    end

    methods (Abstract)        
        create_vars(obj, order, reps); %create the variables of the iqc
        build_psi(obj, vars, order, reps); %build the iqc
        build_M(obj, vars, reps); %build the running cost
        build_X(obj, vars, order, reps); %build the terminal cost
        build_loop(obj, reps); %form the signal transformation matrix

        filter_constraints(obj, cons, order, vars, rho_sched, iqc); %form the constraints      that must be satsified by the filter coefficients
    end

end