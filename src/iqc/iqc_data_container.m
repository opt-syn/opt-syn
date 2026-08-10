classdef iqc_data_container
    %IQC_DATA container for the IQC arising from operators
    
    properties
        iqc; % (iqc_loop_split/iqc_loop_factored)the iqc
        iqc_op; %accumulated iqc for operators
        m_same; %(bool) for Sml, does m=L? If so, there is no uncertainty
        ind_same; %(int array) indices for which m = L
        task; %(string) analysis or synthesis
        augmented=false; %(bool) true if reduced-order control;
        ERGODIC = false; %(bool) true if function values
        rotate = false; %(bool) true if periodic-orbit
        common_rho = 1; %(double) common rho used in noncausal synthesis
    end
    
    methods
        function obj = iqc_data_container()
            %IQC_DATA Constructor            
        end
        
    end
end

