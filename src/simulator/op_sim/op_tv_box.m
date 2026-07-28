classdef op_tv_box < op_sim_interface
    %OP_TV_BOX total variation and box constraint
        

    properties
        BOX=1; %size of box [0, BOX]
        lam_tv; %regularization for total variation
        m; %horizontal dimension
        n; %vertical dimension
    end
    
    methods
        function obj = op_tv_box(BOX, lam_tv, m, n)
            %OP_TV_BOX constructor for l1 norm ball
            %   operations used in the evaluation of the operator
                        
          
            
            obj@op_sim_interface();

            if nargin 
                obj.BOX = BOX;
                obj.lam_tv = lam_tv;
                obj.m = m;
                obj.n = n;
            end
            
        end          

         function w = fw(obj, k, z, param)
            %forward evaluation of the procedure oracle w = E'(E z- b) 
            %
            %Args: 
            %   k (int): time index
            %   z:       input to oracle
            %   param:   parameter structure for the operator
            %
            %Returns:
            %   w:       the w such that w = F(z)

            m = obj.m;
            n = obj.n;

            X = reshape(z, m, n, 3);
            G = zeros(m, n, 3);
        
            for c = 1:3
                Xc = X(:,:,c);
        
                [gx, gy] = gradient(Xc);
        
                denom = sqrt(gx.^2 + gy.^2);
                denom = max(denom, 1e-12);   % avoid divide by zero
        
                px = gx ./ denom;
                py = gy ./ denom;
        
                G(:,:,c) = -obj.lam_tv * divergence(px, py);
            end
        
            w = G(:);
        end

        function z = bw(obj, k, D, v, param)
            %backwards evaluation of an oracle, generalization of a 
            %proximal evaluation with preconditioner D                      
            %
            %Args: 
            %   k (int): time index
            %   D:       prox parameter
            %   v:       input to proximal oracle
            %   param:   parameter structure for the operator
            %
            %Returns:
            %   z:       the z such that z = (I - D F)^(-1)(v)


            %coordinate dimensions
            dl = obj.blocksize(v);
            m = obj.m;
            n = obj.n;

            %effective regularization
            alpha = obj.lam_tv * D;

            Z = reshape(v, m, n, 3);

            tau = 0.25;
            n_inner = 50;
            % Dskron = kron(eye(dl), sqrt(D));
            
            X = zeros(m, n, 3);
            for c = 1:3
                px = zeros(m,n);
                py = zeros(m,n);
                for k = 1:n_inner
                    div_p = divergence(px, py);
                    U = div_p - Z(:,:,c)/alpha;
        
                    [gx, gy] = gradient(U);
        
                    px = px + tau * gx;
                    py = py + tau * gy;
        
                    normp = max(1, sqrt(px.^2 + py.^2));
                    px = px ./ normp;
                    py = py ./ normp;
                end
                X(:,:,c) = Z(:,:,c) - alpha * divergence(px, py);
                X(:,:,c) = min(max(X(:,:,c), 0), 1);
            end
        
            z = X(:);

            
        end


        function f_out = f(obj, k, z, param)
            %total variation penalty. Ignores the box constraint in
            %function evaluation
            %
            %Args: 
            %   k (int): time index
            %   z:       input to oracle           
            %   param:   parameter structure for the operator
            %
            %Returns:
            %   f_out:   f_out = TV(z)

            
            X = reshape(z, obj.m, obj.n, 3);
        
            gx = zeros(obj.m, obj.n, 3);
            gy = zeros(obj.m, obj.n, 3);
            
            for c = 1:3
                gx(1:end-1,:,c) = X(2:end,:,c) - X(1:end-1,:,c);
                gy(:,1:end-1,c) = X(:,2:end,c) - X(:,1:end-1,c);
            end
        
            % Isotropic TV
            tv = sum( sqrt( gx.^2 + gy.^2 ), 'all' );
        
            f_out = obj.lam_tv * tv;
        end
    end
end

