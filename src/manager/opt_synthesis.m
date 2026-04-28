classdef opt_synthesis < opt_manager_interface
    %OPT_SYNTHESIS synthesis of optimization algorithms
    %
    % iterative procedure to find a point beta satisfying
    % the fixed-point equation 
    %               0 \in sum_i F_i(\beta).
    %
    % in which the oracles F_i are interfaced over a dynamical network
    properties
        opts = struct('reduced_order',  false);   
        iqc_rob = {};
    end
    
    methods
        function obj = opt_synthesis(sys)
            %OPT_SYNTHESIS Construct an instance of this class
            %   Detailed explanation goes here

            obj@opt_manager_interface(sys);
            
        end
        
        function obj = process_argument(obj,iqc_rob)
            %PROCESS_ARGUMENT assign orders to the operators/IQCs
            
            %iqc_rob: IQCs representing the robust uncertainties

            if nargin > 1
                obj.iqc_rob = iqc_rob;
            end
        end

        %% extract the solution                   
        function  sol = process_recovery(obj, sol, lmi_out);
            %PROCESS_RECOVERY recover the controller from the solution
            
            error('TODO: Controller Recovery to be implemented')
            
            % iqc_rec = cell(size(obj.iqc_op));
            % for i = 1:length(obj.iqc_op)
            %     if isnumeric(obj.iqc_op{i})
            %         %the Same oracle (m=L, known linear transformation)
            %         iqc_rec{i} = obj.iqc_op{i};
            %     else
            %         iqc_rec{i} = obj.iqc_op{i}.recover(lmi_out);
            %     end
            % 
            % end
            % 
            % sol.iqc = iqc_rec;
        end
    end
end

