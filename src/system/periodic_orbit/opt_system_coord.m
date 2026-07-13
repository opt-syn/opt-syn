classdef  opt_system_coord < opt_system_periodic_orbit
    %OPT_SYSTEM_COORD interconnection of network and operators
    %
    %specifically of a coordinate descent type
    %
    %a periodic system: repeated and predictable cycle evaluation  
    %
    % orbit: the periodicity is highly structured in a symmetric manner
    %
    %
    %a periodic system: repeated and predictable cycle evaluation   
    %
    % w(k) \in F(z(k))
    %
    % [x(k+1)] = [A(k)    Bw(k)    Bwp(k)   Bu(k)  ][x(k)]   state transition
    % [z(k)  ] = [Cz(k)   Dzw(k)   Dzwp(k)  Dzu(k) ][w(k)]   output to oracle
    % [zp(k) ] = [Czp(k)  Dzpw(k)  Dzpwp(k) Dzpu(k)][wp(k)]  output to performance
    % [y(k)  ] = [Cy(k)   Dyw(k)   Dywp(k)  Dyu(k) ][u(k)]   output to controller
    %
    %
    % related together by orbit matrices R describing the symmetry
    %
    % Example:
    %  A(k)  = [Rx 0]^-k [A(0)  Bw(0)]  [Rx 0]^k
    %  Cz(k) = [0 Rz]    [Cz(k) Dzw(0)] [0 Rw]
    

    properties
        % ncoord;
        coord_ind_u;
    end
    
    methods
        function obj = opt_system_coord(op, P, K, coord_ind_u, bind, tracking)
            %OPT_SYSTEM_COORD coordinate   

            %
            %The system P should be in original form.
            %
            %The system K should already be lifted into the coordinate
            %descent form?
            if nargin < 5
                s = length(op);
                bind = 1:s;            
            end

            if nargin < 6
                 tracking = [];
            end

            %which coordinates should be stored a la coordinate descent?            
            cc = length(coord_ind_u);


            %isolate the coordinate descent structure
            Pprim_all = coordinate_descent_primitives(c);
            Pprim = Pprim_all{1};

            P_use = [];
            


            %attach the original system to a coordinate descent subsystem
            Pl = P.lift(cc);            
            Pcoord = lft(Pl, Pprim);
            
            


            obj@opt_system_periodic_orbit(op, P, K, bind, tracking, R);

            
            obj.coord_ind_u = coord_ind_u;
        end        


        function mode_next = next_mode(obj, mode)
            %next mode in the switching sequence
            nc = obj.ncoord;
            mode_next = 1+ mod(mode, nc);
        end

    end
end

