classdef lmi_analysis_periodic < lmi_analysis_interface
    %LMI_ANALYSIS_PERIODIC analysis LMIs for algorithmic interconnections
    %involving periodic linear networks and controllers
    %
    % w(k) \in F(z(k))
    %
    % [x(k+1)] = [A(k)    B(k)     Bp(k)   ][x(k)]   state transition
    % [z(k)  ] = [Cz(k)   Dzw(k)   Dzwp(k) ][w(k)]   output to oracle
    % [zp(k) ] = [Czp(k)  Dzpw(k)  Dzpwp(k)][wp(k)]  output to performance
    %
    %A(k) = A(k+T) for some known time T
    %
    %instances of these algorithms include cyclic coordinate descent
    %methods. Periodic systems can also be unrolled into an LTI system
    %(monodromy methods)


    properties
        opts = struct("COMMON", false);
    end

    methods
        function obj = lmi_analysis_periodic(sys)
            %LMI_DISPATCH_LTI Construct an instance of this class
            %   Detailed explanation goes here
            obj@lmi_analysis_interface(sys);
        end

        function ns = Nss(obj)
            %NSS: Number of subsystems            
            ns = obj.sys.Nss;
        end

        function [vars_diss, cons]= create_vars_storage(obj, alg_psi, cons, name)
            %create_vars_storage create variables for the dissipation
            %constraints

            if nargin < 4
                name = [];
            end

            G_cell = cell(obj.Nss, 1);

            if obj.opts.COMMON
                %common storage function among all subsystems
                [vars_diss, cons] = create_vars_storage@lmi_analysis_interface(obj, alg_psi, cons, name);

                G = vars_diss.G;
                G_cell = cell(obj.Nss, 1);
                for i = 1:obj.Nss
                    G_cell{i} = G;
                end
                
            else
                %define a storage function for each subsystem

                for i = 1:obj.Nss
                    G_curr = obj.define_storage_G(num2str(i));
                    G_cell{i} = G_curr;
                end

            end

        end
    end
end