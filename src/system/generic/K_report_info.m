classdef K_report_info
    %K_REPORT_INFO information about the recovered controller in synthesis
    
    properties
        K_sub; %subcontroller, before internal model
        K;     %controller, after internal model
        model; %the internal model
        alg_trans; %the plant with confirmed performance by LMIs
    end
    
    methods
        function obj = K_report_info(inputArg1,inputArg2)
            %K_REPORT_INFO Constructor                       
        end

    end
end

