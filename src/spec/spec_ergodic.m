classdef spec_ergodic < spec_interface
    %SPEC_ERGODIC enforce asymptotic stability of the algorithmic
    %interconnection (ergodic/sublinear rates)
    
    properties
        type = 'ergodic';
        rho = 1;
        Nw = 1; %weighted consensus matrix
    end
    
    methods
        function [obj, sys_erg]= spec_ergodic(sys)
            %SPEC_ERGODIC Construct an instance of this class
            %   Detailed explanation goes here

            

            c = sys.op{1}.c;
            s = length(sys.op);

            Nw0 = sys.get_consensus_weighted();
            Nw = kron(Nw0, eye(c));


            sys_erg = sys;
            % iwp = [];
            % izp = [];
            % Nw = [];

            %TODO: figure out the duality gap relation in context with
            %multi-oracles
            [P_erg, iwp, izp] = sys.P.perf_ergodic(Nw);

            sys_erg.P = P_erg;
            % sys_erg.nwp = P_erg.nwp;
            % sys_erg.nzp = P_erg.nzp;
            

            obj@spec_interface(iwp, izp);
            
            obj.Nw = Nw;
            
        end

        function M = supply(obj, vars_spec)
            %SUPPLY generate the supply rate
            %
            %duality gap 
            if isempty(obj.Nw)
                M = [];
            else
                % [na, nb] = size(obj.Nw);
                
                % M0 = [0, -1; -1, 0];

                %this carries the duality gap to the other side

                %sum_i f(z_i) - <u*_i, y_i - y*_i> <= Supply

                M0 = [0, 1; 1, 0];
                % M0  = zeros(2);
                Mk = kron(M0, eye(length(obj.izp)));

                outer = blkdiag(eye(length(obj.izp)), obj.Nw);

                M = outer' * Mk * outer;
            end

        end

    end
end

