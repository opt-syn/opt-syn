classdef iqc_loop_interface
    %IQC_INTERFACE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Psi1 = 1;%primal filter (output of nonlinearity)
        Psi2 = 1; %dual filter   (input of nonlinearity);        
        M=0;
        X=0;
        loop = [];
    end
    
    methods
        function obj = iqc_loop_interface(Psi1, M, loop, Psi2, X)
            %IQC_LOOP_INTERFACE Construct an instance of this class
            %   Detailed explanation goes here
            % obj.Property1 = inputArg1 + inputArg2;

            obj.Psi1 = Psi1;
            obj.M = M;

            if nargin > 2
                obj.loop = loop;
            else
                np = size(obj.Psi1.D, 1);
                nq = size(M, 1) - np;
                obj.loop = [zeros(np, np), eye(nq,np);
                        eye(np, nq), zeros(nq, nq)];
            end

            if nargin > 3
                obj.Psi2 = Psi2;
                obj.X = X;
            end
        end
        
        % function outputArg = method1(obj,inputArg)
        %     %METHOD1 Summary of this method goes here
        %     %   Detailed explanation goes here
        %     outputArg = obj.Property1 + inputArg;
        % end
    end

    methods (Abstract)
        nf(obj)
        get_psi(obj)
    end
end

