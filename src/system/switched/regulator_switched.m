classdef regulator_switched < regulator_interface
    %REGULATOR_SWITCHED Regulator for switched systems 
    %
    % [x(k+1)] = [A(mode(k))    Bd(mode(k))    Bu(mode(k))  ][x(k)]   state transition
    % [e(k)  ] = [Ce(mode(k))   Ded(mode(k))   Deu(mode(k)) ][d(k)]   output to  regulated error    
    % [u(k) ] =  [Cy(mode(k))   Dyd(mode(k))   Dyu(mode(k)) ][u(k)]   output to controller    


    %The system obeys a switching logic: restricted mode transitions 
    %from mode(k) to mode(k+1) based on a switching graph (adjacency matrix
    % sys.adj)
    %
    %instances of these algorithms include optimization algorithms under
    %a-priori-unknown time-varying delays

    methods
        function obj = regulator_switched(sys)
            %REGULATOR_SWITCHED build the regulator
            
            obj@regulator_interface(sys)
        end

        function ns = Nss(obj)
            %NSS: number of states
            ns = obj.sys.Nss;
        end

        function obj = form_internal_model(obj)
            %FORM_INTERNAL_MODEL create the internal model by solving the regulator
            %equation. Inputs are the system (P, bind, tracking, op)
            %
            %op is important for which oracles are equaltiy constarints and
            %which are inequality constraints

            Npre = obj.sys.get_consensus(obj.sys.op, obj.sys.bind);
            c = obj.sys.op{1}.c; %coordinate lifts: change this later?
            N = kron(Npre, eye(c));

            [sN0, dN0] = size(Npre);
            [sN, dN] = size(N);
            n = obj.sys.P.nx;
            
            nu = obj.sys.P.nu;

            [src, dst] = obj.sys.get_arcs();
            Narcs = length(src);

            [Sbeta, Rbeta] = obj.sys.get_tracked_opt();

            if isempty(obj.sys.tracking)
                S = eye(size(N, 2)+c);
                R = S;
                reg_ans = [];
                reg_mat = [];

                
                ns = size(S, 1);
                %go through each subsystem
                for i = 1:Narcs
                    Pcurr = obj.sys.P{src(i)};
                    [A, B1, B2, C1, D11, D12, C2, D21, D22] = Pcurr.ss_zy_wu();

                    reg_ans_curr = [zeros(n, c), -B1*N;  -kron(ones(sN0,1), eye(c)), -D11*N];
                    reg_mat_curr = [A, B2; C1, D12];


                    reg_mat_next = zeros(size(reg_mat_curr));
                    reg_mat_next(1:(c*n), 1:(c*n)) = -eye(c*n);

                    sz2 = size(reg_mat_curr, 2);
                    reg_mat_embed = zeros(size(reg_mat_curr, 1), sz2*obj.Nss);

                    shift_curr = (1:sz2) + (src(i)-1)*sz2;
                    
                    reg_mat_embed(:, shift_curr) = reg_mat_curr;

                    if src(i) == dst(i)
                        reg_mat_embed(:, shift_curr) = reg_mat_embed(:, shift_curr) + reg_mat_next;
                    else
                        shift_next = (1:sz2) + (dst(i)-1)*sz2;
                        reg_mat_embed(:, shift_next) = reg_mat_next;
                    end

                    % reg_mat_embed(:, shift_curr) = reg_mat_curr + reg_mat_next;

                    reg_ans = [reg_ans; reg_ans_curr];

                    
                    reg_mat = [reg_mat; reg_mat_embed];
                end
              
                %solve the regulator equation
                null_basis = null(reg_mat, 'rational');
                % try                    
                    sol0 = lsqminnorm(reg_mat, reg_ans);
                % catch
                sol_err = reg_mat * sol0 - reg_ans;
                if norm(sol_err) > 1e-12
                    error('Regulator equation cannot be solved')
                end

                nnull = size(null_basis, 2);

                %extract the solution

                Pi0 = cell(obj.Nss, 1);
                Gam0 = cell(obj.Nss, 1);
                Phi0 = cell(obj.Nss, 1);

                Pi_basis = cell(obj.Nss, 1);
                Gam_basis = cell(obj.Nss, 1);
                Phi_basis = cell(obj.Nss, 1);

                count = 0;
                for i = 1:obj.Nss
                    %get the regulator equation solution
                    Pcurr = obj.sys.P{i};
                    [~, ~, ~, ~, ~, ~, C2, D21, D22] = Pcurr.ss_zy_wu();


                    ind_pi = count + (1:n);
                    ind_gam = count + n+ (1:nu);
                    Pi0{i} = sol0(ind_pi, :);
                    Gam0{i} = sol0(ind_gam, :);
                    Phi0{i} = D21 * [zeros(sN, c), N] + D22*Gam0{i} + C2*Pi0{i};

                    %get the free parameters
                    if nnull
                        Pi_basis_pre = null_basis(ind_pi, :);
                        Gam_basis_pre = null_basis(ind_gam, :);
                        Phi_basis_pre = D22*Gam_basis_pre + C2*Pi_basis_pre;

                        Pi_basis{i} = kron(Pi_basis_pre, eye(nnull));
                        Gam_basis{i} = kron(Gam_basis_pre, eye(nnull));
                        Phi_basis{i} = kron(Phi_basis_pre, eye(nnull));
                    else
                        Pi_basis{i} = [];
                        Gam_basis{i} = [];
                        Phi_basis{i} = [];
                    end

                    count = count + n + nu;

                end                                    
                

            else
                error('Switched regulation: tracking not yet supported')
            end


            obj.S = S;
            obj.R = R;
            obj.Pi = Pi0;
            obj.Gam = Gam0;
            obj.Phi = Phi0;
            obj.Pi_basis = Pi_basis;
            obj.Gam_basis = Gam_basis;
            obj.Phi_basis = Phi_basis;
        end

        function [regulator_closed] = check_regulator(obj)
            %CHECK_REGULATOR is the regulator equation satisfied?
            sys_cl = lft(obj.sys.P, obj.sys.K);

            % error('Switched regulator: regulator check not yet done')
            Npre = obj.sys.get_consensus(obj.sys.op, obj.sys.bind);
            c = size(sys_cl{1}.D, 1)/length(obj.sys.bind); %coordinate lifts: change this later?
            N = kron(Npre, eye(c));


            [sN0, dN0] = size(Npre);
            [sN, dN] = size(N);


            [src, dst] = obj.sys.get_arcs();
            Narcs = length(src);

            if isempty(obj.sys.tracking)
                S = eye(size(N, 2)+1);
                R = S;
                reg_ans = [];
                reg_mat = [];

                n = length(sys_cl{1}.A);

                %go through each subsystem


                for i = 1:Narcs
                    Pcurr = obj.sys.P{src(i)};
                    [A, B, C, D] = ssdata(sys_cl.P{src(i)});


                    reg_ans_curr = [zeros(n, 1), -B*N;  -ones(sN, 1), -D*N];
                    reg_mat_curr = [A; C];
                    reg_mat_next = [-eye(n); zeros(size(C))];

                    sz2 = size(reg_mat_curr, 2);
                   
                    reg_mat_embed = zeros(size(reg_mat_curr, 1), sz2*obj.Nss);

                    shift_curr = (1:sz2) + (src(i)-1)*sz2;

                    reg_mat_embed(:, shift_curr) = reg_mat_curr;

                    if src(i) == dst(i)
                        reg_mat_embed(:, shift_curr) = reg_mat_embed(:, shift_curr) + reg_mat_next;
                    else
                        shift_next = (1:sz2) + (dst(i)-1)*sz2;
                        reg_mat_embed(:, shift_next) = reg_mat_next;
                    end

                    % reg_mat_embed(:, shift_curr) = reg_mat_curr + reg_mat_next;

                    reg_ans = [reg_ans; reg_ans_curr];


                    reg_mat = [reg_mat; reg_mat_embed];
                end

                %solve the regulator equation
                null_basis = null(reg_mat, 'rational');
                
                sol0 = reg_mat \ reg_ans;
                sol_err = reg_mat * sol0 - reg_ans;
                if norm(sol_err) > 1e-12
                    error('Regulator equation cannot be solved')
                end
            
                nnull = size(null_basis, 2);

                %extract the solution

                Pi0 = cell(obj.Nss, 1);
                Th0 = cell(obj.Nss, 1);               

                count = 0;
                nxn = obj.sys.nxn;                
                nxi = obj.sys.nxi;
                for i = 1:obj.Nss
                    %get the regulator equation solution                    
                    ind_pi = count + (1:nxn);
                    ind_th = count + nxn + (1:nxi);
                    Pi0{i} = sol0(ind_pi, :);
                    Th0{i} = sol0(ind_th, :);                    

                    count = count + nxn + nxi;

                end                                    
            end

            regulator_closed = struct('S', S, 'R', R, 'Pi', Pi0, ...
                'Th', Th0 );
        

        end

        %% use the model in synthesis

        function sys = get_model(obj, ind, vars_reg)
            %get_model
            %fetch the internal model (nominal) at mode 'ind'
            %
            %
            %with edits: allow for selection of model within feasible set

            %TODO: allow for parameterizations based on the variables

            
            if nargin < 3
                Phi = obj.Phi{ind};
                Gam = obj.Gam{ind};
            else
                Phi = vars_reg.Phi{ind};
                Gam = vars_reg.Gam{ind};
            end

            [nu, ns] = ssize(Gam);
            ny = ssize(Phi, 1);


            Am = obj.S;
            Bm = [zeros(ns, ny), eye(ns), zeros(ns, nu)];
            Cm = [-Gam; Phi];
            Dm = [zeros(nu, ny), zeros(nu, ns), eye(nu);
                eye(ny), zeros(ny, ns), zeros(ny, nu)];

            n = struct;
            n.nw = ny;
            n.nz = nu;
            n.nu = nu + ns;
            n.ny = ny;
            n.nzp = 0;
            n.nwp = 0;

            P = ss(Am, Bm, Cm, Dm, 1);
            sys = genplant(P, n);

        end

        function plant_model = connect_model(obj, plant, ind, rho)
            %connect the model (nominal regulator equation)

            if nargin < 4
                rho = 1;
            end

            if iscell(plant)
                plant_model = cell(obj.Nss, 1);
                for i = 1:obj.Nss
                    model = obj.get_model(i);
                    model_rho = rhotrafo(model, rho);
                    plant_model{i} = lft(plant{i}, model_rho);
                end
            else
                model = obj.get_model(ind);
                model_rho = rhotrafo(model, rho);
                plant_model = lft(plant, model_rho);
            end

            
        end



    end
end