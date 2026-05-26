classdef spec_ergodic < spec_interface
    %SPEC_ERGODIC enforce asymptotic stability of the algorithmic
    %interconnection (ergodic/sublinear rates)
    
    properties
        type = 'ergodic';
        rho = 1;
        erg_pass = 0;
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

        function [quad, objective] = supply_quad(obj, vars_spec)
            %SUPPLY_QUAD decomposed quadratic performance specification

            [quad, objective] = supply_quad@spec_interface(obj, vars_spec);

            if obj.target
                %add a strict passivity penalty 
                nwp = length(obj.iwp);
                Q_new = drep(vars_spec.erg_pass, nwp);

                quad.Q = quad.Q + Q_new;
                % quad.Q = quad.Q - Q_new;

                objective  = vars_spec.erg_pass;
                
            end
        end


        function [obj] = set_p(obj, p)
            %SET_P set a parameter when performing bisection
            %
            %
            %Example: Peak-to-Peak norm certifier
            %or l2 Gain bound

            obj.erg_pass = p;

        end


        function [vars, cons] = create_vars(obj, cons, name, config)
            %CREATE_VARS form the variables for the problem                        
            if nargin < 3
                name = [];
            end

            if config.LMILAB
                erg_pass = lmim(['erg_pass', name], 1, 1, 'sym');            
            else
                erg_pass = sdpvar(1, 1);            
            end

            vars = struct('erg_pass', erg_pass);
            if obj.target
                cons = append_lmi(cons, -erg_pass+ 100, config.LMILAB);
                cons = append_lmi(cons, erg_pass+ 100, config.LMILAB);
            end
        end
    

        function M = supply(obj, vars_spec)
            %SUPPLY generate the supply rate
            %
            %duality gap 
            if isempty(obj.Nw)
                M = [];
            else
                %passivity relation (kind of)
                M0 = [0, 1; 1, 0];                
                Mk = kron(M0, eye(length(obj.izp)));

                outer = blkdiag(eye(length(obj.izp)), obj.Nw);

                M = outer' * Mk * outer;


                %add `strict' passivity penalty, this offers an objective
                %in alternating methods for ergodic convergence
                nzp = length(obj.izp);
                nwp = length(obj.iwp);

                M_strict = blkdiag(zeros(nzp), eye(nwp)*obj.erg_pass);

                M = M + M_strict;

            end

        end

    end
end

