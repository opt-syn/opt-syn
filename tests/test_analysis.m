classdef test_analysis < matlab.unittest.TestCase

    methods(TestClassSetup)
        % Shared setup for the entire test class
        function addFolders(~)
            addpath(genpath('../src'));
        end
    end

    methods(TestMethodSetup)
        % Setup for each test
    end

    methods(Test)
        % Test methods

        function testAnalysisSML(testCase)
            %analysis of gradient-descent algorithm with optimal step size
            m= 1;
            L = 10;

            %different operator classes for op1
            op1_sml = op_sml(m, L);
            op1_quad = op_quad(m, L);

            % gradient descent
            gamma = 2/(L+m);
            sK = ss([1], [-gamma], [1], [0],1);

            sys_sml = opt_system({op1_sml}, [],  sK);
            sys_quad = opt_system({op1_quad}, [],  sK);

            order = [2];

            man_sml = opt_analysis(sys_sml);
            man_quad = opt_analysis(sys_quad);

            sol_sml = man_sml.bisect(order);
            rho_sml = sol_sml.rho;

            sol_quad = man_quad.bisect(order);
            rho_quad = sol_quad.rho;


            optimal_rate = (L-m)/(L+m); % 0.8182
            %% analyze solution
            %what is this?
            % P_sml = op1_sml.dhd_lift(order, sol_sml.vars.op{1}, sol_sml.cert.iqc_op{1});


            % should we also check regulator equation?

            testCase.verifyEqual(round(rho_sml, 3), round(optimal_rate, 3));
            testCase.verifyEqual(round(rho_quad, 3), round(optimal_rate, 3));

        end
    end

end