classdef test_simulation < matlab.unittest.TestCase

    methods(TestClassSetup)
        function addFolders(~)
            addpath(genpath('../src'));
        end
    end

    methods(TestMethodSetup)
        % Setup for each test
    end

    methods(Test)
        % Test methods

        function testPGDsim(testCase)
            rng(32, 'twister');

            %PGD to minimize quadratic under hard l1 ball constraint

            d = 400; %dimension of variable beta

            %define the quadratic
            m = 1; L = 1000;
            Q = rand_quad(d, m, L);
            zstar = 100*(2*rand(d, 1) - 1);
            op1 = op_sim_quad(Q, zstar);


            %define the L1 ball
            tau = 200;
            op2 = op_sim_l1_hard(tau);
            ops = {op1, op2};

            %PGD algorithm
            % gamma = 2/(L + m);
            gamma = 1/L;
            K = ss(1, [-gamma, -gamma], [1; 1], [0, 0; -gamma, -gamma],1);

            %form the system
            sys = opt_system(ops, [], K);

            %simulate and plot
            sim = alg_sim(sys, d);
            T = 20;
            sim.sampler.x0 = 10*randn(1, d);
            sim_out= sim.sim(T);
            res_w = sim_out.res_w;
            res_z = sim_out.res_z;
            f = sim_out.f;

            testCase.verifyTrue(res_w(end)< res_w(1))
            testCase.verifyTrue(res_z(end) < res_z(1))
            testCase.verifyTrue(f(end) < f(1))
        end
    end

end