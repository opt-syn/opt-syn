classdef op_sml_variational < op_sml_interface
    properties
        rho_val; % Decay rate / discount factor used in the filter
    end
    
    methods
        function obj = op_sml_variational(m, L, c)
            if nargin < 4
                c = 1;
            end
            obj@op_sml_interface(m, L, c);
        end
        
        function [vars] = create_vars(obj, order, reps)
            % sector multiplier: lambda_sec
            % off-by-one multiplier: lambda_off
            
            % note: no theory for order >1, so order=1 always
            % reps: I still dont get how the interface defines repititions,
            % so I leave reps=1 for now

            lambda_sec = lmim(['lambda_sec_', obj.sid], 1, 1, 'full');
            lambda_off = lmim(['lambda_off_', obj.sid], 1, 1, 'full');
            vars = struct('lambda_sec', lambda_sec, 'lambda_off', lambda_off);
        end
        
        function [Psi1, Psi2] = build_psi(obj, vars, order, reps)

            m_val = obj.m;
            L_val = obj.L;
            a = sqrt(m_val * (L_val - m_val) / 2);
            
            % State-space matrices 
            A_block = zeros(4);
            B_block = [1,  0,  1,  0;
                       0,  1,  0, -1;
                       a,  0,  0,  0;
                       -m_val, 1,  0,  0];
                       
            C_block = [-L_val, 1, 0, 0;
                       0, 0, 0, 0;
                       0, 0, 1, 0;
                       a, 0, 0, 0;
                       0, 0, 0, 1;
                       -m_val, 1, 0, 0];
                       
            D_block = [L_val, -1, 0, 0;
                       -m_val, 1, 0, 0;
                       zeros(4, 4)];

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Input 3 to the filter is x_t^* - x_{t+1}^* = C (xi_t^* - xi_{t+1}^*) 
            % where x are the algo iterates and xi the algo state. So the
            % correct filter matrices should be
            % A_block = zeros(4);
            % B_block = [1,  0,      C,             0;
            %            0,  1,      zeros(1,nxi), -1;
            %            a,  0,      zeros(1,nxi),  0;
            %            -m_val, 1,  zeros(1,nxi),  0];
            % 
            % C_block = [-L_val, 1, 0, 0;
            %            0, 0, 0, 0;
            %            0, 0, 1, 0;
            %            a, 0, 0, 0;
            %            0, 0, 0, 1;
            %            -m_val, 1, 0, 0];
            % 
            % D_block = [L_val, -1, zeros(1,nxi), 0;
            %            -m_val, 1, zeros(1,nxi), 0;
            %            zeros(4, 3+nxi)];
            % 
            % Since the output matrix C depends on the synth controller, I 
            % dont know how to access C here, and how to convexify it for
            % the synthesis. I therefore keep x as input to the third
            % channel and hope that its possible later to left/right 
            % multiply a blkdiag(I, I, C, I) somewhere.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                       
            % Generalize to multiple repetitions using block diagonal
            Af = kron(eye(reps), A_block);
            Bf = kron(eye(reps), B_block);
            Cf = kron(eye(reps), C_block);
            Df = kron(eye(reps), D_block);
            
            Psi1 = sdpss(Af, Bf, Cf, Df);
            Psi2 = ss(eye(reps)); % no acausal stuff here
        end
        
        function M = build_M(obj, vars, order, reps)
            M2 = [0, 1; 1, 0];
            M6 = blkdiag(0.5*[0, 1; 1, 0], [1, 0; 0, -1], 0.5*[1, 0; 0, -1]);

            sec_val = lmim_index(vars.lambda_sec, 1, 1);
            off_val = lmim_index(vars.lambda_off, 1, 1);
            
            M_sec = sec_val * M2;
            M_off = off_val * M6;
            
            M = blkdiag(M_sec, M_off);
        end
        
        function X_out = build_X(obj, vars, order, reps)
            X_out = []; % no terminal constraint
        end
        
        function loop = build_loop(obj, reps)
            % var IQC is implemeted without loop transformation
            loop = eye(2 * reps);
        end
        
        function cons = filter_constraints(obj, cons, order, vars, rho_sched, iqc)
            % nonegativity of lambdas is only constraint
            cons = elem_nonneg(vars.lambda_sec, cons, obj.LMILAB);
            cons = elem_nonneg(vars.lambda_off, cons, obj.LMILAB);
        end
    end
end