classdef alg_sim
    %ALG_SIM_INTERFACE execution of the algorithmic interconnection
    
    properties
        sys; %system to simulate
        d;   %number of dimensions (kronecker lift)        
        c=1;   %number of partitions of dimension
        sampler = struct('wp', @(param) [], 'param', @(param) [], ...
            'x0', @() [], 'param0', @() []);
    end
    
    methods 
        function obj = alg_sim(sys, d, c, sampler)
            %ALG_SIM_INTERFACE Construct an instance of this class
            %   Detailed explanation goes here
            obj.sys = sys;
            obj.d = d;
            
            if nargin >=3
                obj.c = c;
            end

            if nargin >= 4
                obj.sampler = sampler;
            end
        end
        
        function ssim = sim(obj, T)
            %SIM: simulate a trajectory execution

            %process the input

            
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

            param0 = obj.sampler.param0;
      

            
            
            %log the signals
            ssim = struct;
            ssim.w = zeros( obj.sys.P.nw, dl, T);
            ssim.wp = zeros( obj.sys.P.nwp, dl, T);
            ssim.z = zeros( obj.sys.P.nz, dl, T);
            ssim.zp = zeros( obj.sys.P.nzp, dl, T);
            ssim.u = zeros( obj.sys.P.nu, dl, T);
            ssim.y = zeros( obj.sys.P.ny, dl, T);
            ssim.xn = zeros( obj.sys.nxn, dl, T);            
            ssim.xi = zeros( obj.sys.nxi, dl, T);         
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
                    wp = obj.sampler.wp(param);
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
                        %use backward evaluation
                        vi_vec = reshape(vi, [], 1);
                        zi_vec = op_curr.bw(k, -Dzw_curr, vi_vec, param);
                        zi = reshape(zi_vec, [], dl);
                        wi = -(Dzw_curr) \ (vi - zi);
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
                %TODO: enable this (and debug it)
                [y, u] = obj.sys.get_internal_signals(param, x, [w; wp]);

                %log the signals


                xn = x(1:obj.sys.nxn, :);
                xi = x((obj.sys.nxn+1):end, :);
                ssim.xn(:, :, k) = xn;
                ssim.xi(:, :, k) = xi;
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
            
            %identify equality constraints
            EQ_mask = logical(cellfun(@(q) q.EQUALITY, obj.sys.op));
            EQ_sum = sum(EQ_mask);         
            nop = length(EQ_mask);
            % sel_z = kron(sparse(1:(nop-EQ_sum), find(~EQ_mask), ones((nop-EQ_sum), 1), (nop-EQ_sum), nop), speye(c));

            sel_w_EQ_pre = kron(sparse(1:EQ_sum, find(EQ_mask), ones(EQ_sum, 1), EQ_sum, nop), speye(c));
            sel_w_EQ = full(sel_w_EQ_pre(:, obj.sys.bind));

            %get the residuals
            %optimality residual: sum(w) = 0            
            wsum2 = pagemtimes(repmat(kron(ones(1, s), eye(c)), 1, 1, T),  ssim.w);                       
            ssim.res_w = sqrt(squeeze(sum(wsum2.^2, [1, 2])))';

            % p_res = pagemtimes(full(sel_w_EQ), ssim.w);
            if EQ_sum 
                w_eq = pagemtimes(repmat(kron(sel_w_EQ, eye(c)), 1, 1, T),  ssim.w);                       
                ssim.eq = sqrt(squeeze(sum(w_eq.^2, [2])))';
            else
                ssim.eq = [];
            end

            %consensus residual: z - zavg = 0           
            % z_sel = pagemtimes(full(sel_z),  ssim.z);
            zavg = pagemtimes(repmat(kron(ones(1, s), eye(c)/s), 1, 1, T),  ssim.z);
            
            ssim.res_z = sqrt(squeeze(sum((ssim.z - repmat(zavg, s, 1, 1)).^2, [1, 2])))';

            ssim.k = 0:(T-1);
        end
    end

end



