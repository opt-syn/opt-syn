classdef cert_synthesis
    %CERT_SYNTHESIS certificate of the synthesis program
    %   undefined

    properties
        K_sub; %subcontroller (before the internal model interconnection)
        K; %controller (after internal model interconnection)
        model; % (genplant/genplant_poly) internal model used to certify regulation constraint
        
        diss; % (diss_data) dissipation constraints
        iqc_op; %(cell of iqc) the IQCs for the individual operators
        iqc_op_all; %(iqc_data) block-diagonal concatenation of the IQCs
        alg_psi; %(genplant/genplant_poly) closed-loop system that is LMI-verified to satisfy performance specifications        
        
        gain; %validation of performance criteria (passivity, H infinity)
        Gcl; %closed-loop storage matrix for alg_psi
        Ycl; %similarity transformation for convexification
        
    end

    methods
        function obj = cert_synthesis()
            %CERT_SYNTHESIS constructor
            
        end
    end
end