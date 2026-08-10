classdef lmi_analysis_periodic_orbit < lmi_analysis_lti
    %LMI_ANALYSIS_PERIODIC_ORBIT analysis LMIs for algorithmic interconnections
    %involving periodic linear networks and controllers
    %
    %Note: this calls routines from lmi_analysis_lti with no extra
    %functionality
    
    % w(k) \in F(z(k))
    %
    % [x(k+1)] = [A(k)    B(k)     Bp(k)   ][x(k)]   state transition
    % [z(k)  ] = [Cz(k)   Dzw(k)   Dzwp(k) ][w(k)]   output to oracle
    % [zp(k) ] = [Czp(k)  Dzpw(k)  Dzpwp(k)][wp(k)]  output to performance
    %
    %A(k) = A(k+T) for some known time T
    %furthermore, matrices [R] are known with 
    %A(k) = (R^k)' A(0) R^k,  R^h = R (and the same for other channels.
    %
    % this is a specialization of general periodic algorithms
    %
    %
    %
    %instances of these algorithms include cyclic coordinate descent
    %methods. Periodic-orbit systems can also be unrolled into an LTI system
    %(monodromy methods): a single large LMI system rather than multiple 
    % coupled smaller LMI systems
    


    properties
        M; %periodicity in the state in dynamics

    end

    methods
        function obj = lmi_analysis_periodic_orbit(sys, config)
            %LMI_ANALYSIS_PERIODIC_ORBIT Construct or.
            %   use the LTI analysis code w.r.t. the periodic-rotated
            %   system and IQCs for the time-varying nonlinearities
            obj@lmi_analysis_lti(sys, config);
            obj.M = sys.M;
        end
 

    end
end