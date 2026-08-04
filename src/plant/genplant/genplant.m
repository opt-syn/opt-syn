classdef genplant
    % GENPLANT A generalized plant with structured channel partitioning.
    %
    % A generalized plant wraps a state-space system ``P`` and partitions
    % its input/output channels into three groups:
    %
    % - **Outputs (from network):** :math:`[z,\; z_p,\; y]` — operator inputs,
    %   performance outputs, and controller inputs.
    % - **Inputs (to network):** :math:`[w,\; w_p,\; u]` — operator outputs,
    %   performance inputs, and controller outputs.
    %
    % The channel dimensions are stored as properties and used by the
    % indexing methods to extract submatrices of ``P``.
    %
    % .. note::
    %
    %    The ordering convention is ``[z, zp, y]`` for outputs and
    %    ``[w, wp, u]`` inputs of the plant.
    %
    % Example::
    %
    %    n.nz = 2; n.nw = 2; n.ny = 1; n.nu = 1; n.s = 1;
    %    G = genplant(ss(A, B, C, D, 1), n);

    properties
        P       % State-space system (``ss`` or ``sdpss`` object).

        s  = 0  % Number of operators.
        nz = 0  % Dimension of operator input :math:`z` (from network).
        nzp = 0 % Dimension of performance output :math:`z_p` (from network).
        ny = 0  % Dimension of controller input :math:`y` (from network).

        nw = 0  % Dimension of operator output :math:`w` (to network).
        nwp = 0 % Dimension of performance input :math:`w_p` (to network).
        nu = 0  % Dimension of controller output :math:`u` (to network).
    end

    methods
        function obj = genplant(P, n)
            % GENPLANT Construct a generalized plant.
            %
            % Wraps a state-space system and, optionally, records the
            % dimensions of the channel partition.
            %
            % :param P: State-space system to wrap.
            % :type P: ss or sdpss
            % :param n: Channel dimensions. Required fields: ``nz``, ``nw``,
            %    ``ny``, ``nu``. Optional fields: ``nzp``, ``nwp``, ``s``.
            % :type n: struct
            % :returns: A new ``genplant`` object.
            % :rtype: genplant

            if isa(P, 'tf')
                P = ss(P);
            end
            obj.P = P;
            if isnumeric(P) && ~isempty(P)
                % static system                
                obj.P = ss(P);
            end
            
            if nargin > 1
                if ~isfield(n, 's')
                    obj.s = 0;
                else
                    obj.s = n.s;
                end
                obj.nz = n.nz;
                obj.nw = n.nw;
                obj.ny = n.ny;
                obj.nu = n.nu;

                if isfield(n, 'nzp')
                    obj.nzp = n.nzp;
                end
                if isfield(n, 'nwp')
                    obj.nwp = n.nwp;
                end
            else
                [obj.nz, obj.nw] = size(P.D); 
            end
            
        end

        %% Channel indexers
        % These methods return index vectors into the rows or columns of
        % the system matrices, respecting the ``[z, zp, y]`` / ``[w, wp, u]``
        % ordering.

        function w_ind = index_w(obj)
            % INDEX_W Indices of the operator-output channel :math:`w`.
            %
            % :returns: Index vector ``1:nw``.
            % :rtype: double (row vector)
            w_ind = 1:obj.nw;
        end

        function wp_ind = index_wp(obj)
            % INDEX_WP Indices of the performance-input channel :math:`w_p`.
            %
            % :returns: Index vector for :math:`w_p`, offset by ``nw``.
            % :rtype: double (row vector)
            wp_ind = obj.nw + (1:obj.nwp);
        end

        function u_ind = index_u(obj)
            % INDEX_U Indices of the controller-output channel :math:`u`.
            %
            % :returns: Index vector for :math:`u`, offset by ``nw + nwp``.
            % :rtype: double (row vector)
            u_ind = obj.nw + obj.nwp + (1:obj.nu);
        end

        function y_ind = index_y(obj)
            % INDEX_Y Indices of the controller-input channel :math:`y`.
            %
            % :returns: Index vector for :math:`y`, offset by ``nz + nzp``.
            % :rtype: double (row vector)
            y_ind = obj.nz + obj.nzp + (1:obj.ny);
        end

        function zp_ind = index_zp(obj)
            % INDEX_ZP Indices of the performance-output channel :math:`z_p`.
            %
            % :returns: Index vector for :math:`z_p`, offset by ``nz``.
            % :rtype: double (row vector)
            zp_ind = obj.nz + (1:obj.nzp);
        end

        function z_ind = index_z(obj)
            % INDEX_Z Indices of the operator-input channel :math:`z`.
            %
            % :returns: Index vector ``1:nz``.
            % :rtype: double (row vector)
            z_ind = 1:obj.nz;
        end

        function wr_ind = index_notu(obj)
            % INDEX_NOTU Indices of the combined :math:`[w,\; w_p]` channels.
            %
            % Returns all input indices except the controller output.
            %
            % :returns: Index vector ``1:(nw + nwp)``.
            % :rtype: double (row vector)
            wr_ind = 1:(obj.nw + obj.nwp);
        end

        function T = Ts(obj)
            % TS Sample time of the wrapped system.
            %
            % :returns: Sample time (``0`` for continuous, ``> 0`` for discrete).
            % :rtype: double
            T = obj.P.Ts;
        end

        %% Extract state-space submatrices
        % Methods named ``Bw``, ``Bu``, etc. extract the columns/rows of
        % ``B``, ``C``, ``D`` corresponding to specific channels.

        function B = Bw(obj)
            % BW Input matrix from operator output :math:`w` to state.
            %
            % Extracts the columns of ``B`` indexed by :mat:meth:`index_w`.
            %
            % :returns: Submatrix :math:`B_w`.
            % :rtype: double
            B0 = obj.B;
            B = B0(:, obj.index_w());
        end

        function B = Bwp(obj)
            % BWP Input matrix from performance input :math:`w_p` to state.
            %
            % Extracts the columns of ``B`` indexed by :mat:meth:`index_wp`.
            %
            % :returns: Submatrix :math:`B_{w_p}`.
            % :rtype: double
            B0 = obj.B;
            B = B0(:, obj.index_wp());
        end

        function B = Bu(obj)
            % BU Input matrix from controller output :math:`u` to state.
            %
            % Extracts the columns of ``B`` indexed by :mat:meth:`index_u`.
            %
            % :returns: Submatrix :math:`B_u`.
            % :rtype: double
            B0 = obj.B;
            B = B0(:, obj.index_u());
        end

        function C = Cz(obj)
            % CZ Output matrix from state to operator input :math:`z`.
            %
            % Extracts the rows of ``C`` indexed by :mat:meth:`index_z`.
            %
            % :returns: Submatrix :math:`C_z`.
            % :rtype: double
            C0 = obj.C;
            C = C0(obj.index_z(), :);
        end

        function C = Czp(obj)
            % CZP Output matrix from state to performance output :math:`z_p`.
            %
            % Extracts the rows of ``C`` indexed by :mat:meth:`index_zp`.
            %
            % :returns: Submatrix :math:`C_{z_p}`.
            % :rtype: double
            C0 = obj.C;
            C = C0(obj.index_zp(), :);
        end

        function C = Cy(obj)
            % CY Output matrix from state to controller input :math:`y`.
            %
            % Extracts the rows of ``C`` indexed by :mat:meth:`index_y`.
            %
            % :returns: Submatrix :math:`C_y`.
            % :rtype: double
            C0 = obj.C;
            C = C0(obj.index_y(), :);
        end

        function D = Dzw(obj)
            % DZW Feedthrough from operator output :math:`w` to operator input :math:`z`.
            %
            % :returns: Submatrix :math:`D_{zw}`.
            % :rtype: double
            iw = obj.index_w();
            iz = obj.index_z();
            D0 = obj.P.D;
            D = D0(iz, iw);
        end

        function D = Dzu(obj)
            % DZU Feedthrough from controller output :math:`u` to operator input :math:`z`.
            %
            % :returns: Submatrix :math:`D_{zu}`.
            % :rtype: double
            iu = obj.index_u();
            iz = obj.index_z();
            D0 = obj.P.D;
            D = D0(iz, iu);
        end

        function D = Dyw(obj)
            % DYW Feedthrough from operator output :math:`w` to controller input :math:`y`.
            %
            % :returns: Submatrix :math:`D_{yw}`.
            % :rtype: double
            iw = obj.index_w();
            iy = obj.index_y();
            D0 = obj.P.D;
            D = D0(iy, iw);
        end

        function D = Dywp(obj)
            % DYWP Feedthrough from performance input :math:`w_p` to controller input :math:`y`.
            %
            % :returns: Submatrix :math:`D_{yw_p}`.
            % :rtype: double
            iw = obj.index_wp();
            iy = obj.index_y();
            D0 = obj.P.D;
            D = D0(iy, iw);
        end

        function D = Dzpwp(obj)
            % DZPWP Feedthrough from performance input :math:`w_p` to performance output :math:`z_p`.
            %
            % :returns: Submatrix :math:`D_{z_pw_p}`.
            % :rtype: double
            iw = obj.index_wp();
            izp = obj.index_zp();
            D0 = obj.P.D;
            D = D0(izp, iw);
        end

        function D = Dyu(obj)
            % DYU Direct feedthrough from controller output :math:`u` to controller input :math:`y`.
            %
            % :returns: Submatrix :math:`D_{yu}`.
            % :rtype: double
            iu = obj.index_u();
            iy = obj.index_y();
            D0 = obj.P.D;
            D = D0(iy, iu);
        end

        function D = Dzwp(obj)
            % DZWP Feedthrough from performance input :math:`w_p` to operator input :math:`z`.
            %
            % :returns: Submatrix :math:`D_{zw_p}`.
            % :rtype: double
            iwp = obj.index_wp();
            iz = obj.index_z();
            D0 = obj.P.D;
            D = D0(iz, iwp);
        end

        function D = Dzpu(obj)
            % DZPU Feedthrough from controller output :math:`u` to performance output :math:`z_p`.
            %
            % :returns: Submatrix :math:`D_{z_pu}`.
            % :rtype: double
            iu = obj.index_u();
            izp = obj.index_zp();
            D0 = obj.P.D;
            D = D0(izp, iu);
        end

        function D = Dzpw(obj)
            % DZPW Feedthrough from operator output :math:`w` to performance output :math:`z_p`.
            %
            % :returns: Submatrix :math:`D_{z_pw}`.
            % :rtype: double
            iw = obj.index_w();
            izp = obj.index_zp();
            D0 = obj.P.D;
            D = D0(izp, iw);
        end

        function [A, B, C, D] = ssdata(obj)
            % SSDATA Extract state-space matrices of the wrapped system.
            %
            % :returns: ``[A, B, C, D]`` — the state-space quadruple.
            % :rtype: double (multiple outputs)
            [A, B, C, D] = ssdata(obj.P);
        end

        function Ao = A(obj)
            % A State matrix of the wrapped system.
            %
            % :returns: Matrix :math:`A`.
            % :rtype: double
            Ao = obj.P.A;
        end

        function Bo = B(obj)
            % B Input matrix of the wrapped system.
            %
            % :returns: Matrix :math:`B`.
            % :rtype: double
            Bo = obj.P.B;
        end

        function Co = C(obj)
            % C Output matrix of the wrapped system.
            %
            % :returns: Matrix :math:`C`.
            % :rtype: double
            Co = obj.P.C;
        end

        function Do = D(obj)
            % D Feedthrough matrix of the wrapped system.
            %
            % :returns: Matrix :math:`D`.
            % :rtype: double
            Do = obj.P.D;
        end

        function sys_drop = drop_performance(obj)
            % DROP_PERFORMANCE Remove the performance channel from the plant.
            %
            % Constructs a new ``genplant`` with ``nzp = 0`` and ``nwp = 0``
            % by extracting only the :math:`[z, y]` / :math:`[w, u]`
            % subsystem.
            %
            % :returns: Plant without the performance channel.
            % :rtype: genplant
            [Aa, B1, B2, C1, D11, D12, C2, D21, D22] = ss_zy_wu(obj);

            A = Aa;
            B = [B1, B2];
            C = [C1; C2];
            D = [D11, D12; D21, D22];
            n = obj.dump_dim;
            n.nzp = 0;
            n.nwp = 0;

            sys_drop = genplant(ss(A, B, C, D, obj.Ts), n);
        end

        function [Aa, B1, B2, C1, D11, D12, C2, D21, D22] = ss_zy_wu(obj)
            % SS_ZY_WU Extract plant matrices for the :math:`[w, u] \to [z, y]` subsystem.
            %
            % Returns the nine matrices of the two-port partition:
            %
            % .. math::
            %
            %    \begin{bmatrix} z \\ y \end{bmatrix}
            %    =
            %    \begin{bmatrix} D_{11} & D_{12} \\ D_{21} & D_{22} \end{bmatrix}
            %    \begin{bmatrix} w \\ u \end{bmatrix}
            %    +
            %    \begin{bmatrix} C_1 \\ C_2 \end{bmatrix} x
            %
            % :returns: ``[A, B1, B2, C1, D11, D12, C2, D21, D22]``
            % :rtype: double (multiple outputs)
            Aa = obj.A;
            B1 = obj.Bw;
            B2 = obj.Bu;
            C1 = obj.Cz;
            C2 = obj.Cy;
            D11 = obj.Dzw;
            D12 = obj.Dzu;
            D21 = obj.Dyw;
            D22 = obj.Dyu;
        end

        function nxo = nx(obj)
            % NX Number of states in the wrapped system.
            %
            % :returns: State dimension.
            % :rtype: int
            nxo = length(obj.P.A);
        end

        function P_out = ss(obj)
            % SS Extract the raw state-space object.
            %
            % :returns: The underlying ``ss`` system.
            % :rtype: ss
            P_out = obj.P;
        end

        function P_out = tf(obj)
            % TF Transfer function representation of the plant.
            %
            % :returns: Transfer function of ``P``.
            % :rtype: tf
            P_out = ss2tf(obj.ss());
        end

        function obj = rhotrafo(obj, rho)
            % RHOTRAFO Apply an exponential discount (rho-transformation).
            %
            % Scales ``A`` and ``B`` by :math:`\rho^{-1}`, which corresponds
            % to an exponential weighting of the signals in discrete time.
            %
            % :param rho: Discount factor.
            % :type rho: double
            % :returns: The transformed plant (modified in place).
            % :rtype: genplant
            obj.P.A = obj.P.A * (rho^(-1));
            obj.P.B = obj.P.B * (rho^(-1));
        end

        %% Overloaded interconnection operations

        function b_out = blkdiag(obj, b2)
            % BLKDIAG Block-diagonal interconnection of two generalized plants.
            %
            % Constructs a new ``genplant`` whose system matrix is the
            % block-diagonal of ``obj.P`` and ``b2.P``, with input and
            % output indices interleaved so that the :math:`[z, z_p, y]` /
            % :math:`[w, w_p, u]` channel structure is preserved.
            %
            % :param b2: The second plant.
            % :type b2: genplant
            % :returns: The block-diagonal plant.
            % :rtype: genplant

            b_out = obj;
            b_out.nw = obj.nw + b2.nw;
            b_out.nwp = obj.nwp + b2.nwp;
            b_out.nz = obj.nz + b2.nz;
            b_out.nzp = obj.nzp + b2.nzp;
            b_out.nu = obj.nu + b2.nu;
            b_out.ny = obj.ny + b2.ny;
            b_out.s = obj.s + b2.s;

            P_diag = blkdiag(obj.P, b2.P);
            nin = obj.nw + obj.nwp + obj.nu;
            nout = obj.nz + obj.nzp + obj.ny;

            ind_in = [1:obj.nw, nin + (1:b2.nw), ...
                      obj.nw + (1:obj.nwp), nin + b2.nw + (1:b2.nwp), ...
                      obj.nw + obj.nwp + (1:obj.nu), nin + b2.nw + b2.nwp + (1:b2.nu)];

            ind_out = [1:obj.nz, nout + (1:b2.nz), ...
                      obj.nz + (1:obj.nzp), nout + b2.nz + (1:b2.nzp), ...
                      obj.nz + obj.nzp + (1:obj.ny), nout + b2.nz + b2.nzp + (1:b2.ny)];

            b_out.P = P_diag(ind_out, ind_in);
        end

        function b_out = lft(obj, b2)
            % LFT Lower linear fractional transformation (star product).
            %
            % Closes the feedback loop between ``obj`` (upper block) and
            % ``b2`` (lower block) along the shared :math:`(u, y)` channels.
            % The resulting plant inherits the operator channels of ``obj``
            % and the controller channels of ``b2``.
            %
            % :param b2: The lower plant or controller.
            % :type b2: genplant or ss
            % :returns: The interconnected plant.
            % :rtype: genplant

            b_out = obj;
            if isa(b2, 'genplant')
                b_out.nw = obj.nw;
                b_out.nwp = obj.nwp + b2.nwp;
                b_out.nz = obj.nz;
                b_out.nzp = obj.nzp + b2.nzp;
                b_out.nu = b2.nu;
                b_out.ny = b2.ny;
                b_out.s = obj.s + b2.s;

                b_out.P = lft(obj.P, b2.P, obj.nu, obj.ny);
            else
                [nu2, ny2] = size(b2.D);
                b_out.P = lft(obj.P, b2, nu2, ny2);
                b_out.nu = b_out.nu - nu2;
                b_out.ny = b_out.ny - ny2;
            end
        end

        function b_out = lft_lower(obj, b2)
            % LFT_LOWER Lower linear fractional transformation (alias).
            %
            % Equivalent to :mat:meth:`lft`. Provided for symmetry with
            % :mat:meth:`lft_upper`.
            %
            % :param b2: The lower plant.
            % :type b2: genplant
            % :returns: The interconnected plant.
            % :rtype: genplant
            b_out = obj.lft(obj, b2);
        end

        function b_out = lft_upper(obj, b2, nz2, nw2)
            % LFT_UPPER Upper linear fractional transformation.
            %
            % Closes the feedback loop with ``b2`` as the upper block and
            % ``obj`` as the lower block. The resulting plant inherits the
            % operator channels of ``b2`` and the controller channels of
            % ``obj``.
            %
            % :param b2: The upper plant.
            % :type b2: genplant or ss
            % :param nz2: Number of output channels to close (used when ``b2``
            %    is a plain ``ss``).
            % :type nz2: int
            % :param nw2: Number of input channels to close (used when ``b2``
            %    is a plain ``ss``).
            % :type nw2: int
            % :returns: The interconnected plant.
            % :rtype: genplant

            b_out = obj;
            if isa(b2, 'genplant')
                b_out.nw = b2.nw;
                b_out.nwp = obj.nwp + b2.nwp;
                b_out.nz = b2.nz;
                b_out.nzp = obj.nzp + b2.nzp;
                b_out.nu = obj.nu;
                b_out.ny = obj.ny;
                b_out.s = obj.s + b2.s;

                b_out.P = lft(b2.P, obj.P, obj.nu, obj.ny);
            else
                [nz2_orig, nw2_orig] = size(b2.D);
                b_out.P = lft(b2, obj.P, nz2, nw2);
                b_out.nz = nz2_orig - nz2;
                b_out.nw = nw2_orig - nw2;
            end
        end

        function obj = lift(obj, c)
            % LIFT Kronecker lift of the plant by an identity matrix.
            %
            % Replaces every matrix :math:`M` in the state-space
            % representation with :math:`M \otimes I_c`. All channel
            % dimensions are multiplied by ``c``.
            %
            % :param c: Lift dimension (size of the identity block).
            % :type c: int
            % :returns: The lifted plant.
            % :rtype: genplant
            Ad = kron(obj.P.A, eye(c));
            Bd = kron(obj.P.B, eye(c));
            Cd = kron(obj.P.C, eye(c));
            Dd = kron(obj.P.D, eye(c));

            obj.P = ss(Ad, Bd, Cd, Dd, 1);

            obj.nz  = obj.nz  * c;
            obj.nzp = obj.nzp * c;
            obj.nw  = obj.nw  * c;
            obj.nwp = obj.nwp * c;
            obj.nu  = obj.nu  * c;
            obj.ny  = obj.ny  * c;
        end

        function n = dump_dim(obj)
            % DUMP_DIM Return channel dimensions as a struct.
            %
            % :returns: Struct with fields ``nw``, ``nwp``, ``nu``, ``ny``,
            %    ``nz``, ``nzp``, ``s``.
            % :rtype: struct
            n = struct('nw', obj.nw, 'nwp', obj.nwp, ...
                'nu', obj.nu, 'ny', obj.ny, ...
                'nz', obj.nz, 'nzp', obj.nzp, ...
                's', obj.s);
        end

        %% Performance channel manipulation

        function [obj, iwp] = add_oracle_input(obj, ind_w, ind_z)
            % ADD_ORACLE_INPUT Add external perturbation inputs at the operator.
            %
            % Introduces additional performance inputs that perturb the
            % operator channels:
            %
            % .. math::
            %
            %    w + \delta w \in F(z + \delta z)
            %
            % The new inputs are appended to the :math:`w_p` channel.
            % No extra outputs are added.
            %
            % :param ind_w: Indices within the :math:`w` channel to perturb.
            % :type ind_w: double (vector)
            % :param ind_z: Indices within the :math:`z` channel to perturb.
            % :type ind_z: double (vector)
            % :returns: ``[obj, iwp]`` — updated plant and new performance
            %    input indices.
            % :rtype: genplant, double

            nwpnew = length(ind_w);
            nzpnew = length(ind_z);

            B = obj.B;
            D = obj.D;
            if nwpnew
                Ew = full(sparse(ind_w, 1:nwpnew, ones(nwpnew, 1), obj.nw, nwpnew));
                Bwnew = B(:, obj.index_w) * Ew;
                Dwnew = D(:, obj.index_w) * Ew;
            else
                Bwnew = [];
                Dwnew = [];
            end
            if nzpnew
                Ez = full(sparse(ind_z, 1:nzpnew, ones(nzpnew, 1), obj.nz, nzpnew));
                Bznew = [];
                Dznew = 1 * [Ez; zeros(nzpnew, obj.nz)];
            else
                Bznew = [];
                Dznew = [];
            end

            B_left  = B(:, [obj.index_w, obj.index_wp]);
            B_right = B(:, [obj.index_u]);
            D_left  = D(:, [obj.index_w, obj.index_wp]);
            D_right = D(:, [obj.index_u]);

            Aold = obj.P.A;
            Bnew = [B_left, Bwnew, Bznew, B_right];
            Cold = obj.P.C;
            Dnew = [D_left, Dwnew, Dznew, D_right];

            obj.P = ss(Aold, Bnew, Cold, Dnew, 1);

            iwp = obj.nwp + (1:(nwpnew + nzpnew));
            obj.nwp = obj.nwp + nwpnew + nzpnew;
        end

        function [obj, iwp] = add_oracle_shift(obj, c)
            % ADD_ORACLE_SHIFT Add a shifted perturbation input at the operator.
            %
            % Introduces performance inputs corresponding to a uniform
            % shift across operator channels:
            %
            % .. math::
            %
            %    w \in F(z + \delta z \otimes \mathbf{1}_s)
            %
            % :param c: Coordinate dimension (default: 1).
            % :type c: int
            % :returns: ``[obj, iwp]`` — updated plant and new performance
            %    input indices.
            % :rtype: genplant, double

            if nargin == 1
                c = 1;
            end

            dl = obj.nz / c;
            nzpnew = dl;
            ind_z = 1:dl;

            B = obj.B;
            D = obj.D;

            if nzpnew
                Ez = kron(ones(dl, 1), eye(c));
                Dznew = 1 * [Ez; zeros(nzpnew, c)];
                Bznew = zeros(size(B, 1), size(Dznew, 2));
            else
                Bznew = [];
                Dznew = [];
            end

            B_left  = B(:, [obj.index_w, obj.index_wp]);
            B_right = B(:, [obj.index_u]);
            D_left  = D(:, [obj.index_w, obj.index_wp]);
            D_right = D(:, [obj.index_u]);

            Aold = obj.P.A;
            Bnew = [B_left, Bznew, B_right];
            Cold = obj.P.C;
            Dnew = [D_left, Dznew, D_right];

            obj.P = ss(Aold, Bnew, Cold, Dnew, 1);

            iwp = obj.nwp + (1:c);
            obj.nwp = obj.nwp + c;
        end

        function [obj, izp] = perf_output_w(obj, iw)
            % PERF_OUTPUT_W Add performance outputs tracking the :math:`w` channel.
            %
            % Appends new rows to the :math:`z_p` channel that directly
            % observe selected operator outputs.
            %
            % :param iw: Indices within :math:`w` to observe.
            % :type iw: double (vector)
            % :returns: ``[obj, izp]`` — updated plant and new performance
            %    output indices.
            % :rtype: genplant, double

            A = obj.P.A;
            B = obj.P.B;
            C = obj.P.C;
            D = obj.P.D;

            nnew = length(iw);

            Ctop = C([obj.index_z(), obj.index_zp()], :);
            Dtop = D([obj.index_z(), obj.index_zp()], :);
            Cbot = C([obj.index_y()], :);
            Dbot = D([obj.index_y()], :);

            n = length(A);
            Czp = zeros(nnew, n);
            Ez = full(sparse(iw, 1:nnew, ones(nnew, 1), n, nnew));

            Dzp = Ez;

            Cnew = [Ctop; Czp; Cbot];
            Dnew = [Dtop; Dzp; Dbot];

            obj.P = ss(A, B, Cnew, Dnew, 1);

            izp = obj.nzp + (1:nnew);
            obj.nzp = obj.nzp + nnew;
        end

        function [obj, iwp, izp] = perf_ergodic(obj, Nw)
            % PERF_ERGODIC Add performance channels for ergodic convergence.
            %
            % Appends performance inputs and outputs that encode an
            % ergodic convergence condition using the consensus matrix
            % ``Nw``.
            %
            % :param Nw: Consensus weighting matrix.  Pass ``[]`` to skip.
            % :type Nw: double or empty
            % :returns: ``[obj, iwp, izp]`` — updated plant, new performance
            %    input indices, and new performance output indices.
            % :rtype: genplant, double, double

            if ~isempty(Nw)
                [nzp_erg, nwp_erg] = size(Nw);

                iwp = obj.nwp + (1:nwp_erg);
                izp = obj.nzp + (1:nzp_erg);

                B = obj.B;
                C = obj.C;
                D = obj.D;

                C_new = C(obj.index_z, :);
                B_new_left   = B(:, [obj.index_w, obj.index_wp]);
                B_new_right  = B(:, [obj.index_u]);
                D_new_left   = D(obj.index_z, [obj.index_w, obj.index_wp]);
                D_new_right  = D(obj.index_z, [obj.index_u]);

                B_new_center = zeros(length(obj.A), nwp_erg);
                D_new_center = D(obj.index_z, obj.index_u) * (-Nw);

                B_new = [B_new_left, B_new_center, B_new_right];
                D_new = [D_new_left, D_new_center, D_new_right];

                D_old_topleft  = D([obj.index_z, obj.index_zp], [obj.index_w, obj.index_wp]);
                D_old_topright = D([obj.index_z, obj.index_zp], obj.index_u);
                D_old_botleft  = D(obj.index_y, [obj.index_w, obj.index_wp]);
                D_old_botright = D(obj.index_y, obj.index_u);

                D_old_top = [D_old_topleft, zeros(obj.nz + obj.nzp, nwp_erg), D_old_topright];
                D_old_bot = [D_old_botleft, zeros(obj.ny, nwp_erg), D_old_botright];

                B_merge = B_new;
                C_merge = [C([obj.index_z, obj.index_zp], :); C_new; C(obj.index_u, :)];

                D_merge = [D_old_top; D_new; D_old_bot];

                obj.P = ss(obj.A, B_merge, C_merge, D_merge, 1);

                obj.nzp = obj.nzp + nzp_erg;
                obj.nwp = obj.nwp + nwp_erg;
            else
                iwp = [];
                izp = [];
            end
        end

        function [obj, izp] = perf_output_opt(obj, c, bind)
            % PERF_OUTPUT_OPT Add performance outputs for the optimality condition.
            %
            % Appends a performance output that sums selected operator
            % outputs:
            %
            % .. math::
            %
            %    z_p = \sum_{i=1}^{s} w^i
            %
            % :param c: Coordinate / Kronecker lift dimension (default: 1).
            % :type c: int
            % :param bind: Indices identifying repeated nonlinearity
            %    evaluations (default: ``1:nw/c``).
            % :type bind: double (vector)
            % :returns: ``[obj, izp]`` — updated plant and new performance
            %    output indices.
            % :rtype: genplant, double

            if nargin < 2
                c = 1;
            end
            if nargin < 3
                bind = (1:obj.nw / c);
            end

            [u, uind] = unique(bind);

            A = obj.P.A;
            B = obj.P.B;
            C = obj.P.C;
            D = obj.P.D;

            ind_w = obj.index_w();
            nnew = length(ind_w);

            Ctop = C([obj.index_z(), obj.index_zp()], :);
            Dtop = D([obj.index_z(), obj.index_zp()], :);
            Cbot = C([obj.index_y()], :);
            Dbot = D([obj.index_y()], :);

            n = length(A);
            s = obj.nw / c;
            Czp = zeros(c, n);

            w_ind_keep = double(ismember(1:s, uind));
            Ezp = kron(w_ind_keep, eye(c));
            Dzp = [Ezp, zeros(size(Ezp, 1), obj.nwp + obj.nu)];

            Cnew = [Ctop; Czp; Cbot];
            Dnew = [Dtop; Dzp; Dbot];

            obj.P = ss(A, B, Cnew, Dnew, 1);

            izp = obj.nzp + (1:c);
            obj.nzp = obj.nzp + c;
        end

        function [obj, izp] = perf_output_z(obj, ind_z)
            % PERF_OUTPUT_Z Add performance outputs tracking the :math:`z` channel.
            %
            % Appends new rows to the :math:`z_p` channel that observe
            % selected operator inputs.
            %
            % :param ind_z: Indices within :math:`z` to observe.
            % :type ind_z: double (vector)
            % :returns: ``[obj, izp]`` — updated plant and new performance
            %    output indices.
            % :rtype: genplant, double

            A = obj.P.A;
            B = obj.P.B;
            C = obj.P.C;
            D = obj.P.D;

            nnew = length(ind_z);

            Ctop = C([obj.index_z(), obj.index_zp()], :);
            Dtop = D([obj.index_z(), obj.index_zp()], :);
            Cbot = C([obj.index_y()], :);
            Dbot = D([obj.index_y()], :);

            n = length(A);

            Ez = full(sparse(1:nnew, ind_z, ones(nnew, 1), ...
                length(ind_z), nnew + obj.nzp + obj.ny));

            Czp = Ez * C;
            Dzp = Ez * D;

            Cnew = [Ctop; Czp; Cbot];
            Dnew = [Dtop; Dzp; Dbot];

            obj.P = ss(A, B, Cnew, Dnew, 1);

            izp = obj.nzp + (1:nnew);
            obj.nzp = obj.nzp + nnew;
        end

        function [obj, izp] = perf_output_con(obj, c, iz)
            % PERF_OUTPUT_CON Add performance outputs for consensus tracking.
            %
            % Appends outputs measuring the deviation from the average:
            %
            % .. math::
            %
            %    z_p^i = z^i - \mathrm{average}(z)
            %
            % :param c: Coordinate / Kronecker lift dimension (default: 1).
            % :type c: int
            % :param iz: Indices of operator inputs to track (default:
            %    ``1:nz``).
            % :type iz: double (vector)
            % :returns: ``[obj, izp]`` — updated plant and new performance
            %    output indices.
            % :rtype: genplant, double

            if nargin < 2
                c = 1;
            end
            if nargin < 3
                iz = 1:obj.nz;
            end

            A = obj.P.A;
            B = obj.P.B;
            C = obj.P.C;
            D = obj.P.D;

            nnew = length(iz);

            Ctop = C([obj.index_z(), obj.index_zp()], :);
            Dtop = D([obj.index_z(), obj.index_zp()], :);
            Cbot = C([obj.index_y()], :);
            Dbot = D([obj.index_y()], :);

            n = length(A);

            Ez = full(sparse(1:nnew, iz, ones(nnew, 1), ...
                length(iz), nnew + obj.nzp + obj.ny));

            Iz = eye(obj.nz);
            Jz = ones(obj.nz, obj.nz) / (obj.nz / c);

            Resz = blkdiag((Iz - Jz), eye(obj.nzp + obj.ny));

            Czp = (Ez * Resz) * C;
            Dzp = (Ez * Resz) * D;

            Cnew = [Ctop; Czp; Cbot];
            Dnew = [Dtop; Dzp; Dbot];

            obj.P = ss(A, B, Cnew, Dnew, 1);

            izp = obj.nzp + (1:nnew);
            obj.nzp = obj.nzp + nnew;
        end

    end
end