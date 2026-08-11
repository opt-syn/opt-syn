classdef quad_param
    %QUAD_PARAM describe the quadratic performance specification.
    %Also used to define IQCs
    %
    %
    % :math:`\begin{pmatrix} w_k \\ z_k \end{pmatrix}^\top \begin{pmatrix} Q & S \\ S^\top & T^\top U^{-1} T \end{pmatrix} \begin{pmatrix} w_k \\ z_k \end{pmatrix}

    % [w]' [Q      S      ] [w]
    % [z]  [S' T' U^{-1} T] [z]    
    
    properties
        Q = []; % w squared term
        S = []; % cross term
        U = []; % inverse in z squared, for Schur Complements
        T = []; % outer part of z squared term
    end
    
    methods
        function obj = quad_param(Q, S, U, T)
            %QUAD_PARAM Constructor
            %
            %
            %Note:
            %   3 inputs: (Q, S, R), and R will be decomposed into R = T' U^{-1} T
            %   4 inputs: (Q, S, U, T)
           
            if nargin > 0
                obj.Q = Q;

                if nargin == 1
                    nQ = ssize(Q);
                    obj.S = zeros(nQ, 0);
                    obj.T = [];
                    obj.U = [];
                else
                    obj.S = S;
    
                    if nargin == 3
                        % (Q, S, R)
                        [RqV, RqD] = eig(Rq);
                        eRq = diag(RqD);
                        ind_pos = find(abs(eRq) > 1e-12);
    
                        obj.T = RqV(:, ind_pos)';
                        obj.Q = diag(1./eRq(ind_pos));
    
    
                    elseif nargin == 4
                        % (Q, S, U, T)
                        obj.U = U;
                        obj.T = T;
                    end   
                    
                end
            end
           
        end
        
        function nw_out = nw(obj)
            %number of inputs
            
            nw_out = ssize(obj.Q, 1);
        end

        function nz_out = nz(obj)
            %number of outputs
            
            nz_out = ssize(obj.U, 1);
        end

        function quad_new = blkdiag(obj,other)
            %BLKDIAG block-diagonals of supply rates, respecting partitions
            %
            %Args:
            %   other (quad_param): the other quadratic supply rate
            %   
            %Return:
            %   quad_new (quad_param): the new quadratic supply rate


            
            Q_new = blkdiag(obj.Q, other.Q);
            S_new = blkdiag(obj.S, other.S);
            U_new = blkdiag(obj.U, other.U);            
            T_new = blkdiag(obj.T, other.T);

            quad_new = quad_param(Q_new, S_new, U_new, T_new);
        end

        function R_out = R(obj)
            % get the double output term

            R_out = obj.T' * (obj.U \ obj.T);


        end
    end
end

