classdef alg_sim
    %ALG_SIM algorithm simulator, execution of the algorithmic interconnection
    %
    %Plots the procedure to solve :math:`0 \in \sum_{i=0}^s F_i(\beta^*)`
    
    properties
        sys; %system to simulate
        d;   %number of dimensions (kronecker lift)        
        c=1;  %number of partitions of dimension
        sampler = alg_sim_sampler(); %random sample routines
        EQUALITY = 0; %Is an op_sim_equality object present? True if so.
    end
    
    methods 
        function obj = alg_sim(sys, d, sampler)
            %ALG_SIM Construct an alg_sim object    
            %
            % Args:
            %   sys: system to simulate
            %   d: number of dimensions (multiplicity of kronecker lift)
            %   c: number of partitions of the dimension/coordinate blocks
            %   sampler: random sampler code used in the algorithm execution. 


            % Notes:
            %  the sampler structure is
            %        wp, performance input
            %        x0, initial state
            %        param, parameters at each time
            %        param0, initial parameters
            
            
            
            obj.sys = sys;
            obj.d = d;
            
            
            obj.c = obj.sys.op{1}.c;
            

            if nargin >= 3               
                obj.sampler = sampler;
            else
                obj.sampler = struct('wp', @(k, param) [], 'param', @(k, param) [], ...
            'x0', @() [], 'param0', @() []);
                
                if ismember(obj.sys.type, ['switched', 'periodic', 'periodic_orbit'])
                    %to be improved
                    obj.sampler.param0 = struct('mode', randi([1, obj.sys.Nss], 1));
                    obj.sampler.param = ...
                        @(param_in)struct('mode',  obj.sys.next_mode(param_in.mode));
                end
                
            end

            

            obj.EQUALITY = any(logical(cellfun(@(q) q.EQUALITY, obj.sys.op)));
        end
        
        function ssim = sim(obj, T)
            %SIM: simulate a trajectory execution
            %
            % Args:
            %   T: Number of steps in execution
            %
            % Returns:
            %   ssim (alg_sim_out): struct that stores the algorithm execution. Each output is indexed by time

            

            
            %dimensions            
            s = length(obj.sys.bind);            
            d = obj.d;
            c = obj.c;
            dl = floor(d/c);
            n=obj.sys.nxi;            

            xz = zeros(obj.sys.nxn + obj.sys.nxi, dl);
         
            if isnumeric(obj.sampler.x0)
                if isempty(obj.sampler.x0)
                    x0 = xz;
                else
                    x0 = obj.sampler.x0;
                end
            else
                x0 = obj.sampler.x0();
                if isempty(x0)
                    x0 = xz;
                end
            end

            param0 = obj.sampler.param0();
                 
            
            %log the signals
            ssim = alg_sim_out();
            ssim.s = length(obj.sys.bind);
            ssim.w = zeros( obj.sys.P.nw, dl, T);
            ssim.wp = zeros( obj.sys.P.nwp, dl, T);
            ssim.z = zeros( obj.sys.P.nz, dl, T);
            ssim.zp = zeros( obj.sys.P.nzp, dl, T);
            ssim.u = zeros( obj.sys.P.nu, dl, T);
            ssim.y = zeros( obj.sys.P.ny, dl, T);
            ssim.xn = zeros( obj.sys.nxn, dl, T);            
            ssim.xc = zeros( obj.sys.nxi, dl, T);         
            ssim.mode = zeros(1, T);
            ssim.param = cell(1, T);
           

            %function values            
            nfi = zeros(s, 1);
            z0 = zeros(d, 1);
            for i = 1:s                
                op_curr = obj.sys.get_op(i);
                fi = op_curr.f(0, z0, param0);
                nfi(i) = length(fi);
            end
            nf = sum(nfi);
            ssim.f = zeros(nf, T);


            %iterated signals
            w = zeros(obj.sys.P.nw, dl);
            wp = zeros(obj.sys.P.nwp, dl);
            z = zeros(obj.sys.P.nz, dl);
            zp = zeros(obj.sys.P.nzp, dl);
            u = zeros(obj.sys.P.nu, dl);
            y = zeros(obj.sys.P.ny, dl);
            param = param0;
            x = x0;
            f = zeros(nf, 1);

            %main loop
            for k = 1:T
                %perform current iteration
                alg = obj.sys.get_alg(param);
                [A, B, C, D] = ssdata(alg);
                Cz = C(obj.sys.P.index_z(), :);
                Czp = C(obj.sys.P.index_zp(), :);

                Dzw = D(obj.sys.P.index_z(), obj.sys.P.index_w());
                Dzpwp = D(obj.sys.P.index_zp(), obj.sys.P.index_wp());
                Dzwp = D(obj.sys.P.index_z(), obj.sys.P.index_wp());
                Dzpw = D(obj.sys.P.index_zp(), obj.sys.P.index_w());

                %TODO: allow for performance
                for i = 1:s
                    %ASSERT D is lower-triangular
                    wp = obj.sampler.wp(k, param);
                    if i==1
                        %flush the w and z values
                        w = NaN * w;
                        
                        f = [];
                        z = NaN * z;
                        
                    end
                    i_ind = c*(i-1) + (1:c);
                    i_rem_ind = 1:(c*(i-1));
                    vi = Cz(i_ind, :) * x ;
                    if obj.sys.P.nwp
                        vi = vi + Dzwp(i_ind, :) * wp;
                    end
                    
                    if i > 1
                        vwi = Dzw(i_ind, i_rem_ind) * w(i_rem_ind, :);
                        vi = vi + vwi;
                    end
                    
                    Dzw_curr = Dzw(i_ind, i_ind);
                    op_curr = obj.sys.get_op(i);
                    if any(Dzw_curr,"all")
                        %use backward evaluation (I+ (-D) F)^-1
                        vi_vec = reshape(vi, [], 1);
                        zi_vec = op_curr.bw(k,  -Dzw_curr, vi_vec, param);
                        zi = reshape(zi_vec, [], dl);
                        wi = (Dzw_curr) \ (zi - vi);
                    else
                        %use forward evaluation
                        zi = vi;
                        
                        vi_vec = reshape(vi, [], 1);   
                        zi_vec = vi_vec;
                        wi_vec = op_curr.fw(k, vi_vec, param);
                        wi = reshape(wi_vec, [], dl);
                    end

                    %function evaluation
                    

                    fi = op_curr.f(k, zi_vec, param);



                    %accumulate the w, z, and f vectors
                    w(i_ind, :) = wi;
                    z(i_ind, :) = zi;
                    
                    f_ind = sum(nfi(1:(i-1))) + (1:nfi(i));
                    f(f_ind, :) = fi;
                end

                %get performance outputs
                if obj.sys.P.nzp
                    zp = Czp * x + Dzpw * w;
                    if obj.sys.P.nwp
                        zp = zp + Dzpwp * wp;
                    end
                end

                %extract internal signals
                
                [y, u] = obj.sys.get_internal_signals(param, x, [w; wp]);

                %log the signals


                xn = x(1:obj.sys.nxn, :);
                xc = x((obj.sys.nxn+1):end, :);
                ssim.xn(:, :, k) = xn;
                ssim.xc(:, :, k) = xc;
                ssim.z(:, :, k) = z;
                ssim.zp(:, :, k) = zp;                
                ssim.w(:, :, k) = w;
                ssim.wp(:, :, k) = wp;
                ssim.u(:, :, k) = u;
                ssim.y(:, :, k) = y;
                ssim.f(:, k) = f;
                if obj.sys.Nss==1
                    ssim.mode(k) = 1;
                else
                    ssim.mode(k) = param.mode;
                end
                
                ssim.param{k} = param;
                
                
                %prepare for new iteration
                xnext = A * x + B * [w; wp];                
                
                %TODO: state-dependent transitions?
                paramnext = obj.sampler.param(param);


                x = xnext;
                param = paramnext;

            end

            %% postprocess the run
            
            %identify equality constraints
            EQ_mask = logical(cellfun(@(q) q.EQUALITY, obj.sys.op));
            EQ_sum = sum(EQ_mask);         
            nop = length(EQ_mask);
            % sel_z = kron(sparse(1:(nop-EQ_sum), find(~EQ_mask), ones((nop-EQ_sum), 1), (nop-EQ_sum), nop), speye(c));



            %only one group of equality constraints is permitted
            %get the residuals
            %optimality residual: sum(w) = 0            
            %or in equality constraints: sum(w) \in ran(E'):
            if EQ_sum 
                %TODO: multistep equality constraint methods
                %indices with equality constraints
                sel_EQ_pre = kron(sparse(1:EQ_sum, find(EQ_mask), ones(EQ_sum, 1), EQ_sum, nop), speye(c));
                sel_EQ = full(sel_EQ_pre(:, obj.sys.bind));

                %indices without equality constraints
                sel_std_pre = kron(sparse(1:(nop-EQ_sum), find(~EQ_mask), ones((nop-EQ_sum), 1), (nop-EQ_sum), nop), speye(c));
                sel_std = full(sel_std_pre(:, obj.sys.bind));

                %find the equality constraints
                w_std = pagemtimes(repmat(kron(sel_std, eye(c)), 1, 1, T),  ssim.w);  
                w_perm = squeeze(reshape(permute(w_std, [3, 1, 2]), T, [], 1))';
                op_eq = find(EQ_mask);

                %evaluate the optimality condition sum(w) \in ran(E')
                E = obj.sys.op{op_eq}.E;
                b = obj.sys.op{op_eq}.b;
                lam_rec = E' \ w_perm;
                eq_err = w_perm - E' * lam_rec;
                ssim.res_w = sqrt(squeeze(sum(eq_err.^2, [1])))';

                %evaluate the primal feasibility condition Ez = b
                z_EQ = pagemtimes(repmat(kron(sel_EQ, eye(c)), 1, 1, T),  ssim.w);  
                z_EQ_perm = squeeze(reshape(permute(z_EQ, [3, 1, 2]), T, [], 1))';
                
                eq_res = E * z_EQ_perm - b;
                
                ssim.eq = sqrt(squeeze(sum(eq_res.^2, 1)))';
            else
                ssim.eq = [];
                %TODO: sums over binds
                bind = reshape(obj.sys.bind, [], 1);
                [gg, gc] = groupcounts(bind);
                bind_weight = reshape(1./gg(bind), 1, []);


                wsum2 = pagemtimes(repmat(kron(bind_weight, eye(c)), 1, 1, T),  ssim.w);                       
                ssim.res_w = sqrt(squeeze(sum(wsum2.^2, [1, 2])))';
            end

            %consensus residual: z - zavg = 0           
            %same in all oracles, regardless of equalities
            %TODO: generalize to different consensus constraints
            zavg = pagemtimes(repmat(kron(ones(1, s), eye(c)/s), 1, 1, T),  ssim.z);
            
            ssim.res_z = sqrt(squeeze(sum((ssim.z - repmat(zavg, s, 1, 1)).^2, [1, 2])))';

            ssim.k = 0:(T-1);
        end
    end




end



