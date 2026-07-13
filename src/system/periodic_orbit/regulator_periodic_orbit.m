classdef regulator_periodic_orbit < regulator_switched
    %REGULATOR_PERIODIC_ORBIT Regulator for periodic systems    
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
            %REGULATOR_PERIODIC undefined
            %   undefined

            sys_per = sys.export_periodic();
            obj@regulator_switched(sys_per)
        end

    end
end