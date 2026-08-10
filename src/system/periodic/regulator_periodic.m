classdef regulator_periodic < regulator_switched
    %REGULATOR_PERIODIC Regulator for periodic systems. This is an alias
    %for the regulator_switched()
    %

    %
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
        function obj = regulator_periodic(sys)
            %REGULATOR_PERIODIC constructor
            obj@regulator_switched(sys)
        end
    end
end