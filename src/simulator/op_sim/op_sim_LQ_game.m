classdef op_sim_LQ_game< op_sim_interface
    %OP_SIM_LQ_GAME pseudogradient of a linear quadratic game
    %
    %Agent payoffs are :math:`J_i(x) = \frac{1}{2}x^\top Q_i x + b_i^\top x
    %+ e_i`. The goal is to find a (generalized variational) Nash
    %equilibrium of the game.
    % 
    %
    %
        

    properties
        Q; %quadratic term (cell)
        b; %linear term (cell)
        e; %constant term (cell)
        n; %partition of agents 

        Q_all; %matrix term in pseudogradient
        b_all; %linear term in pseudogradient
    end
    
    methods
        function obj = op_sim_LQ_game(Q, b, e, n)
            %OP_SIM_LQ_GAME constructor for the pseudogradient
            %   operations used in the evaluation of the operator
            %
            %Args:
            %   Q (cell): each agent's quadratic term 
            %   b (cell): each agent's linear term
            %   c (cell): each agent's constant term
            %   n (int):  number of inputs per agent, sum(n) = length(b)

          
            
            obj@op_sim_interface();

           
            obj.Q = Q;
            obj.b = b;
            obj.e = e;

            if nargin < 3
                n = length(e);
            end

            obj.n = n;

            if isscalar(n)                
                obj.n = repmat(n, length(b)/n);
            end
            
            
            %compute the matrices of the game
            N = length(obj.e);
            Q_all = zeros(length(b));
            b_all = zeros(size(b));
            for i = 1:N
                Q_curr = obj.Q{i};
                b_curr = obj.b{i};

                ind_curr = (1:n(i)) + sum(n(1:i-1));
                ind_other = setdiff(1:(sum(n)), ind_curr);

                %relevant parameters in pseudogradient
                Q_self = Q_curr(ind_curr, ind_curr);
                Q_other = Q_curr(ind_curr, ind_other);
                b_self = b_curr(ind_curr);

                %assemble local parameters into a single matrix
                Q_all(ind_curr, ind_curr) = Q_self;
                Q_all(ind_curr, ind_other) = Q_other;
                b_all(ind_curr) = b_self;
            end
        
            obj.Q_all = Q_all;
            obj.b_all = b_all;



        end          



        function w = fw(obj, k, z, param)
            %forward evaluation of the pseudogradient map
            %
            %Args: 
            %   k (int): time index
            %   z:       input to oracle
            %   param:   parameter structure for the operator
            %
            %Returns:
            %   w:       the w such that w = F(z)

            %pseudogradient
            w = obj.Q_all*z + obj.b_all;
        end

        function z = bw(obj, k, D, v, param)
            %backwards evaluation of an pseudogradient, generalization of a 
            %proximal evaluation with preconditioner D                      
            %
            %Args: 
            %   k (int): time index
            %   D:       prox parameter
            %   v:       input to proximal oracle
            %   param:   parameter structure for the operator
            %
            %Returns:
            %   z:       the z such that z = (I + D F)^(-1)(v)

            

            Dkron = kron(eye(obj.blocksize(v)), D);
            z = (obj.Q_all + Dkron) \ ( Dkron * v - obj.b_all);


        end

        function f_out = f(obj, k, z, param)
            %payoff functions of the game
            %
            %Args: 
            %   k (int): time index
            %   z:       input to oracle           
            %   param:   parameter structure for the operator
            %
            %Returns:
            %   f_out:   f_out = [J_1, J_2, ..., J_n]

            
            N = length(obj.b);
            f_out = zeros(N, 1);
            for i = 1:N

                payoff_curr = 0.5*(z'*obj.Q{i}*z) + ...
                    obj.b{i}'*z + obj.e{i};


                f_out(i) = payoff_curr;
    
            end
        end

        function mu = monotone(obj)
            %montonocity constant of the game
            mu = min(eig(obj.Q_all + obj.Q_all'));
        end

        function L = lipschitz(obj)
            %lipschitz constant  of the game
            L = max(eig(obj.Q_all + obj.Q_all'));
        end

        function L = coco(obj)
            %cocoercivity property of the game
            iQ = inv(obj.Q_all);
            L = min(eig( iQ + iQ'));
        end
    end
end

