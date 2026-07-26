classdef cert_analysis
    %CERT_ANALYSIS certificate of the analysis program

    properties
        iqc_op; %(cell of iqc) the IQCs for the individual operators
        alg_psi; %(genplant/genplant_poly) closed-loop system that is LMI-verified to satisfy performance specifications
        diss; % (diss_data) dissipation constraints
    end

    methods
        function obj = cert_analysis()
            %CERT_ANALYSIS constructor
        end

    end
end