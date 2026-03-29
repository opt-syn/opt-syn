classdef (Abstract) operator_interface
    %OPERATOR_INTERFACE An operator with relation w \in F(z)
    % container for the variables, IQCs, and relations    
    
    properties
        vars; %variables in the analysis/synthesis programs
        id;   %the ID of the operator in the problem (index)
        iqc;  %IQC descring relations satisfied by the operator
        %monotone = true; %can the operator be bounded using monotonicity?
    end
    
    methods
        function obj = operator_interface(id)
            %OPERATOR_INTERFACE Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 2
                obj.id = 0;
            end
        end 
    end

    methods (Abstract)
        %constructor
        create_vars(obj, order, reps)
        create_iqc(obj, order, reps)
        %recovery
%         factor_iqc(obj)        
    end

end