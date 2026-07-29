classdef spec_ergodic < spec_interface
    % SPEC_ERGODIC Ergodic (sublinear-rate) convergence specification.
    %
    % Enforces asymptotic stability of the algorithmic interconnection in an
    % ergodic (averaged) sense, appropriate for algorithms whose last iterate
    % need not converge but whose running averages do. The specification is
    % built on the consensus-weighted performance channel produced by
    % :mat:meth:`plant.genplant.genplant.perf_ergodic`, using the weighted
    % consensus matrix :math:`N_w`.
    %
    % A strict-passivity margin ``erg_pass`` supplies an objective for the
    % alternating methods used to certify ergodic convergence.

    properties

        erg_pass = 0; % Strict-passivity margin used as an objective.
        Nw = 1;       % Weighted consensus matrix :math:`N_w`.
    end
    
    methods
        function [obj, sys_erg]= spec_ergodic(sys)
            % SPEC_ERGODIC Construct an ergodic convergence specification.
            %
            % Reads the consensus structure from the system, appends the
            % ergodic performance channel to the plant via
            % :mat:meth:`plant.genplant.genplant.perf_ergodic`, and returns
            % both the specification and the augmented system.
            %
            % :param sys: The optimization-algorithm system.
            % :returns: ``[obj, sys_erg]`` — the ergodic specification and the
            %    system with the ergodic performance channel added.
            % :rtype: spec_ergodic, (system)

            

            c = sys.op{1}.c;
            s = length(sys.op);

            Nw0 = sys.get_consensus_weighted();
            Nw = kron(Nw0, eye(c));


            sys_erg = sys;

            [P_erg, iwp, izp] = sys.P.perf_ergodic(Nw);

            sys_erg.P = P_erg;
            

            obj@spec_interface(iwp, izp);
            
            obj.Nw = Nw;
            obj.type = 'ergodic';        
            
        end

        function [quad, objective] = supply_quad(obj, vars_spec)
            % SUPPLY_QUAD Decomposed quadratic supply rate.
            %
            % Starts from the base-class decomposition and, when this
            % specification is the optimization target, adds a strict
            % passivity penalty ``erg_pass`` to the :math:`Q` block and uses
            % it as the objective.
            %
            % :param vars_spec: Specification variables (expects ``erg_pass``).
            % :type vars_spec: struct
            % :returns: ``[quad, objective]`` — the decomposed quadratic and the
            %    objective contribution.
            % :rtype: quad_param, double

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
            % SET_P Set the bisection parameter (the strict-passivity margin).
            %
            % :param p: New value of the margin, e.g. for a peak-to-peak
            %    certifier or an :math:`\ell_2`-gain bound.
            % :type p: double
            % :returns: The updated specification.
            % :rtype: spec_ergodic

            obj.erg_pass = p;

        end


        function [vars, cons] = create_vars(obj, cons, name, config)
            % CREATE_VARS Create the strict-passivity margin variable.
            %
            % :param cons: Accumulated LMI constraints.
            % :param name: Name suffix for the created variable.
            % :type name: char
            % :param config: Configuration options.
            % :type config: opt_config
            % :returns: ``[vars, cons]`` — struct with field ``erg_pass`` and the
            %    updated constraint set.
            % :rtype: struct, cell
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
            % SUPPLY Ergodic supply rate (duality gap).
            %
            % Builds a passivity-type supply rate from the consensus matrix
            % :math:`N_w` and adds a strict-passivity penalty scaled by
            % ``erg_pass``. Returns empty if :math:`N_w` is empty.
            %
            % :param vars_spec: Specification variables.
            % :type vars_spec: struct
            % :returns: Quadratic running-cost matrix :math:`M`.
            % :rtype: double

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
