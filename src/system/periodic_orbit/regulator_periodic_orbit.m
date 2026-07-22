classdef regulator_periodic_orbit < regulator_lti
    %REGULATOR_PERIODIC_ORBIT Regulator for periodic-orbit systems.
    %Use routines directly from regulator_lti.
    
    %
    % [x(k+1)] = [A(k)    Bd(k)    Bu(k)  ][x(k)]   state transition
    % [e(k)  ] = [Ce(k)   Ded(k)   Deu(k) ][d(k)]   output to  regulated error    
    % [zp(k) ] = [Cy(k)   Dyd(k)   Dyu(k) ][u(k)]   output to controller    

    %A(k) = A(k+T) for some known time T
    %
    %instances of these algorithms include cyclic coordinate descent
    %methods. Periodic systems can also be unrolled into an LTI system
    %(monodromy methods): a single large LMI system rather than multiple 
    % coupled smaller LMI systems

    methods
        function obj = regulator_periodic_orbit(sys)
            %REGULATOR_PERIODIC_ORBIT constructor

            % sys_per = sys.export_periodic();
            obj@regulator_lti(sys)
        end

        function [A, B1, B2, C1, D11, D12, C2, D21, D22] = ss_zy_wu(obj, param)
            %get plant matrices for the system
            if nargin < 2
                param = [];
            end
            [A, B1, B2, C1, D11, D12, C2, D21, D22] = ss_zy_wu@regulator_lti(obj, param);

            %go to the rotating coordinate frame
            c = size(obj.sys.M, 1);
            n = size(A, 1);
            Rkron = kron(eye(n/c), obj.sys.M);

            A = Rkron * A;
            B1 = Rkron * B1;
            B2 = Rkron * B2;
        end        

        function [S, R] = exosystem(obj, param)
            %get the exosystem at each mode/internal model
            if nargin < 2
                param = [];
            end
            [S, R] = exosystem@regulator_lti(obj, param);

            c = size(obj.sys.M, 1);
            d = size(S, 1);
            Rkron = kron(eye(d/c), obj.sys.M);

            %go to the rotating coordinate frame
            S = Rkron * S;

        end

        function Kcurr = get_K(obj, param)
            %get the controller
            if nargin < 2
                param = [];
            end
            Kcurr = obj.sys.get_K(param);

            %go to the rotating coordinate frame
            c = size(obj.sys.M, 1);
            n = size(Kcurr.A, 1);
            Rkron = kron(eye(n/c), obj.sys.M);

            Kcurr.A = Rkron * Kcurr.A;
            Kcurr.B = Rkron * Kcurr.B;

        end



    end
end