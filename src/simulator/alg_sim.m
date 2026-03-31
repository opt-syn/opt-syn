classdef alg_sim
    %ALG_SIM_INTERFACE execution of the algorithmic interconnection
    
    properties
        sys;
        blocksize;
        sampler;
    end
    
    methods 
        function obj = alg_sim(sys, blocksize, sampler)
            %ALG_SIM_INTERFACE Construct an instance of this class
            %   Detailed explanation goes here
            obj.sys = sys;
            obj.blocksize = blocksize;
            obj.sampler = sampler;
        end
        
        function ssim = sim(obj, T, x0, param0)
            %SIM: simulate a trajectory execution

            %process the input

            
            %dimensions            
            s = length(obj.sys.bind);            
            d = sum(obj.blocksize);
            c = length(obj.blocksize);
            n=obj.sys.n;            

            if nargin >= 3
                if length(x0) == n*d
                    x_vec=x0;
                else
                    x_vec=kron(x0,ones(d,1));
                end
            else
                x_vec=zeros(n*d,1);
            end

            if nargin == 4

            end

            
            
            %log the signals
            ssim = struct;
            ssim.w = zeros( obj.sys.nw, d, T);
            ssim.wp = zeros( obj.sys.nwp, d, T);
            ssim.z = zeros( obj.sys.nz, d, T);
            ssim.zp = zeros( obj.sys.nzp, d, T);
            ssim.u = zeros( obj.sys.nu, d, T);
            ssim.y = zeros( obj.sys.ny, d, T);
            ssim.x = zeros( n, d, T);            
            ssim.param = cell(1, T);
           

            %function values            
            nfi = zeros(s, 1);
            z0 = zeros(d, 1);
            for i = 1:s                
                op_curr = obj.sys.get_op(i);
                fi = op_curr.f(z0, param0);
                nfi(i) = length(i);
            end
            nf = sum(nfi);
            ssim.f = zeros(T, nf);


            %iterated signals
            w = zeros(obj.sys.nw, d);
            wp = zeros(obj.sys.nwp, d);
            z = zeros(obj.sys.nz, d);
            zp = zeros(obj.sys.nzp, d);
            u = zeros(obj.sys.nu, d);
            y = zeros(obj.sys.ny, d);
            x = reshape(x_vec, n, d);
            param = param0;
            f = zeros(nf, 1);

            %main loop
            for k = 1:T
                %perform current iteration
                alg = obj.get_alg(param);
                [A, B, C, D] = ssdata(alg);
                Cz = C(obj.P.index_z(), :);
                Czp = C(obj.P.index_zp(), :);

                Dzw = D(obj.P.index_z(), obj.P.index_w());
                Dzpwp = D(obj.P.index_zp(), obj.P.index_wp());
                Dzwp = D(obj.P.index_z(), obj.P.index_wp());
                Dzpw = D(obj.P.index_zp(), obj.P.index_w());

                %TODO: allow for performance
                for i = 1:s
                    %ASSERT D is lower-triangular
                    wp = obj.sampler.w(param);
                    if i==1
                        %flush the w and z values
                        w = NaN * w;
                        
                        f = [];
                        z = NaN * z;
                        
                    end
                    i_ind = c*(i-1) + (1:c);
                    i_rem_ind = 1:(c*(i-1));
                    vi = Cz(i_ind, :) * x + Dzwp(i_ind, :) * wp;
                    if i > 1
                        vwi = Dzw(i_ind, i_rem_ind) * w(i_rem_ind, :);
                        vi = vi + vwi;
                    end
                    
                    Dzw_curr = Dzw(i_ind, i_ind);
                    op_curr = obj.sys.get_op(i);
                    if any(Dzw_curr,"all")
                        %use backward evaluation
                        
                        zi = op_curr.bw(-Dzw_curr, vi, param);
                        wi = -(Dzw_curr) \ (vi - zi);
                    else
                        %use forward evaluation
                        zi = vi;
                        wi = op_curr.fw(vi, param);
                    end

                    %function evaluation
                    

                    fi = op_curr.f(zi, param);



                    %accumulate the w, z, and f vectors
                    w(i_ind, :) = wi;
                    z(i_ind, :) = zi;
                    
                    f_ind = sum(nfi(1:(i-1))) + (1:nfi(i));
                    f(f_ind, :) = fi;
                end

                %get performance outputs
                zp = Czp * x + Dzpwp * wp + Dzpw * w;

                %extract internal signals
                [y, u] = obj.sys.get_internal_signals(param, x, [w; wp], z);

                %log the signals

                ssim.x(:, :, k) = x;
                ssim.z(:, :, k) = z;
                ssim.zp(:, :, k) = zp;                
                ssim.w(:, :, k) = w;
                ssim.wp(:, :, k) = wp;
                ssim.u(:, :, k) = u;
                ssim.y(:, :, k) = y;
                ssim.param{k} = param;
                
                
                %prepare for new iteration
                xnext = A * x + B * [w; wp];                
                
                %TODO: state-dependent transitions?
                paramnext = obj.sampler.param(param);


                x = xnext;
                param = paramnext;

            end

        end
    end

end


