classdef lmi_dispatch_interface < handle
    %LMI_DISPATCH_INTERFACE analysis and synthesis LMIs for the
    %algorithmic interconnections
    %
    %This contains generic routines common among both analysis and
    %synthesis for every system type
    %  
    properties
        sys;
        config;
    end
    
    methods
        function obj = lmi_dispatch_interface(sys, config)
            %LMI_DISPATCH_INTERFACE Construct an instance of this class
            %   Detailed explanation goes here
            obj.sys = sys;
            obj.config = config;
        end
        

        
        %KYP lemma terms, commonly found matrices
        
        function [cons, objective, vars, con_M] = cons_dynamic(obj, vars, cons, diss)
            %CONS form the dissipation and sign constraints
            %
            %Input:
            %   vars:   variables of the problem        
            %   cons:   accumulated constraints
            %   diss:   structure describing the problem
            %       plant:  system to control
            %       spec:   performance specification           
            %       target: whether the performance measure should be optimized
            %               true:  soft constraint (e.g. Schur complement
            %                                       formulation)
            %               false: hard constraint            
            %
            %Output:
            %   cons:   accumulated constraints
            %   objective:  term to be minimized            
            
            
            %need to look up the right constraint            

            %Upper-levels: iterate over the systems
            [cons, objective, con_M] = obj.con_dynamic_single(vars, cons, diss);

            
                      
        end
    

        function [cons, objective, con_M] = con_dynamic_single(obj,  vars, cons, diss)
            %CON form a single dissipation and sign constraint
            %
            %Input:
            %   vars:   variables of the problem            
            %   diss:   structure describing the problem
            %       plant:  system to control
            %       spec:   performance specification           
            %       target: whether the performance measure should be optimized
            %               true:  soft constraint (e.g. Schur complement
            %                                       formulation)
            %               false: hard constraint
            %   param:  other parameters
            %
            %Output:        
            %   objective:  term to be minimized            
    
    
            %need to look up the right constraint            
    
            
    
            if ismember(diss.spec.type, methods(obj))
                lmi_handle = str2func(diss.spec.type);
                [cons, objective, con_M] = lmi_handle(obj, vars, cons, diss);
            else
                spt = obj.sys.type;
                msg = ['The performance specification ', diss.spec.type, ...
                    ' is not supported for ', spt, 'systems.'];
                error('msg', 'OPT:spec_unsupported');
            end
        end

        %% helper functions
        function M = merge_spec_M(obj, iqc_rob, sp,  vars_spec)
            %MERGE_SPEC_M merge the running cost of the robustness and the 
            %performance specification
            %
            %Input:
            %   iqc_rob: robust IQC 
            %   sp:      performance specification
            %   vars_spec: variables in the specification

            if nargin > 2
                supp =sp.supply(vars_spec);
            else
                supp = sp.supply();
            end
            n1 = iqc_rob.np;
            m1 = iqc_rob.nq;
            n2 = length(sp.izp);
            m2 = length(sp.iwp);
            
            [M] = outer_blkdiag(iqc_rob.M, supp, n1, m1, n2, m2);
            % Mdiag = blkdiag(iqc_op.M, sp.supply);
        end

        function [quad] = quad_objective(obj, M_quad, ind_p, ind_q)
            %QUAD_OBJECTIVE untangle the quadratic objective into a
            %linearizable formulation

            %maybe get rid of this, replace by quad_objective_decomp

            %R = T' U^-1 T, R >0            
            %use eigenvalue arguments here

            Qq = M_quad(ind_q, ind_q);
            Sq = M_quad(ind_q, ind_p);
            Rq = M_quad(ind_p, ind_p);


            [RqV, RqD] = eig(Rq);
            eRq = diag(RqD);
            ind_pos = find(abs(eRq) > 1e-12);

            Tq = RqV(:, ind_pos)';
            Uq = diag(1./eRq(ind_pos));

            quad = struct('Q', Qq, 'S', Sq, 'U', Uq, 'T', Tq);
        end

        function quad_m = merge_quad(obj, M_rob, M_spec)
            %merge together quadratic performance specifications;


            Qq= blkdiag(M_rob.Q, M_spec.Q);
            Sq= blkdiag(M_rob.S, M_spec.S);
            Tq= blkdiag(M_rob.T, M_spec.T);
            Uq= blkdiag(M_rob.U, M_spec.U);

            quad_m = struct('Q', Qq, 'S', Sq, 'U', Uq, 'T', Tq);
        end



        function sb = sys_block(obj, plant, Pnew, Pold)
            % SYS_BLOCK system block used in analysis programs
            %
            %sb =  [0, I]^T [Pold, 0] [0, I]
            %      [A, B]   [0,      -Pnew] [A, B]

            
            A = plant.A;
            B = plant.B;
            
            [n, m] = size(B);  

            Ablock = [eye(n), zeros(n, m);
                A, B];

            Pblock = blkdiag(Pold, -Pnew);

            sb = Ablock' * Pblock * Ablock;            

        end

        function sb = supply_block(obj, plant, M)
            % SUPPLY_BLOCK supply block used in analysis programs
            %
            %sb =  [C, D]^T [-M] [C D]
            %               


            Cblock = [plant.C, plant.D];   
            
            sb = Cblock' * (-M) * Cblock;

        end

        function [ind_q, ind_p] = get_idx_performance(obj, diss);
            %GET_IDX_PERFORMANCE 
        
            nz = diss.iqc_rob.nz + diss.spec.nzp;
            nw = diss.iqc_rob.nw + diss.spec.nwp;

            ind_p = 1:nz;
            ind_q = nz + (1:nw);

        end


        function verdict = LMILAB(obj)
            %is LMILAB used?
            verdict = obj.config.LMILAB();
        end


        function [plant_no_p, CDp] = separate_performance_output(obj, diss)

            %SEPARATE_PERFORMANCE_OUTPUT
            %extract the performance output from the plant
            %
            %plant_no_p:    the plant with the performance output removed
            %CDp:           the entries of [Cp, Dp] matrix for the
            %               performance output
            %
            %
            %This routine is used in the computation of l2 norms via Schur
            %complements            
            
            %separate the performance outputs   
            nzp = length(diss.spec.izp);
            ind_sep = (diss.iqc_rob.np) + (1:nzp);
            nz = ssize(diss.plant.D, 1);

            ind_diff_sep = setdiff(1:nz, ind_sep);
            Iz = eye(nz);
            Izp = Iz(ind_diff_sep, :);

            %the plant without the schur-complemented-out performance input
            plant_no_p = Izp*diss.plant;


            Ezp = sparse(1:length(ind_sep), ind_sep, ones(length(ind_sep)), ...
                length(ind_sep), ssize(diss.plant.C, 1));


            CDp = Ezp * [diss.plant.C, diss.plant.D];

        end

        function cons = con_spread(obj, cons, vars)
            %CON_SPREAD increase numerical conditioning by separating the 
            %primal and dual blocks
            % in synthesis only

            %do nothing in analysis           

        end

        function [cons, objective, con_M] = stability(obj, vars, cons, diss)
            %STABILITY certification of exponential stability
            %
            %the supply function in the specification is empty,
            %so just call quadratic performance.


            [cons, objective, con_M] = obj.quad(vars, cons, diss);

        end

        function [cons, objective, con_M] = ergodic(obj, vars, cons, diss)
            %ERGODIC certification of ergodic convergence
            %            
            %call quadratic performance.
            %
            %ergodic means the system should have previously been adjusted


            [cons, objective, con_M] = obj.quad(vars, cons, diss);

        end


        function [cons, objective, con_M] = e2e(obj, vars, cons, diss)
            %E2E: energy to energy gain

            if diss.spec.target
                [cons, objective, con_M] = obj.e2e_target(vars, cons, diss);
            else
                %is a special case of quadratic performance
                [cons, objective, con_M] = obj.quad(vars, cons, diss);
            end           
        end


        
    end




    methods (Abstract)
        create_vars(obj, vars, cons, alg, specs)
        create_vars_storage(obj, cons, alg_psi, name)
        % process_recovery(obj, sol, lmi_out, alg_psi)

        %supported performance measures
        % quad 
        % stability(obj, vars, cons, diss)  
        % e2e(obj, vars, cons, diss)  
        % e2e_target(obj, vars, cons, diss)  

        %should allow for:
            %quad
            %e2e
            %e2p
        % stability(obj, vars, cons, diss)
        
        
    end
end


