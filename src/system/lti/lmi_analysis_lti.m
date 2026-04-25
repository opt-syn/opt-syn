classdef lmi_analysis_lti < lmi_analysis_interface
    %LMI_ANALYSIS_LTI analysis LMIs for algorithmic interconnections
    %involving linear-time-invariant (LTI) networks and controllers
    %
    % w(k) \in F(z(k))
    %
    % [x(k+1)] = [A    B     Bp   ][x(k)]   state transition
    % [z(k)  ] = [Cz   Dzw   Dzwp ][w(k)]   output to oracle
    % [zp(k) ] = [Czp  Dzpw  Dzpwp][wp(k)]  output to performance
    %
    % performance specification: wp -> zp from (spec)
    %
    %   Implemented
    %
    %
    %   TODO:
    %       stability
    %       h2
    %       e2p
    %       p2p
    %
    
    methods
        function obj = lmi_analysis_lti(sys)
            %LMI_DISPATCH_LTI Construct an instance of this class
            %   Detailed explanation goes here
            obj@lmi_analysis_interface(sys);
        end

        %function [vars_diss, cons]= create_vars_storage(obj, alg_psi, cons, name)
        %is the default
        
        
            %CON_TERMINAL terminal constraint
        

        function [vars, cons] = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

