classdef test_synthesis < matlab.unittest.TestCase

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

        function testSynthesis(testCase)
            %describe the operators
            m = 1; L = 50;
            op1 = op_sml(m, L); %gradient of f
            op2 = op_pcc();     %indicator function of L1 ball


            %run the synthesis routine, use bisection to minimize the convergence rate
            sys = opt_system({op1, op2});

            man = opt_synthesis(sys);
            sol = man.bisect();

            rho = sol.rho; % 0.8676
            
            testCase.verifyTrue(rho<1);
            %% with time-delay dynamics
            delay = [1, 0];
            network = bridge_channel_delay(delay, delay);
            sys = opt_system({op1, op2}, network);
            man = opt_synthesis(sys);
            sol_delay = man.bisect();  % 0.9860

            testCase.verifyTrue(sol_delay.rho < 1);
        end
    end

end