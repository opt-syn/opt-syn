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
           
            if nargin > 0
                obj.Q = Q;
                obj.S = S;
                obj.U = U;
                obj.T = T;
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
    end
end

