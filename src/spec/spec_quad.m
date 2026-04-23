classdef spec_quad < spec_interface
    %SPEC_QUAD quadratic performance
    %
    %
    %
    %[zp]' [Q  S]  [zp] > 0
    %[wp]  [S' R]  [wp] 

    
    properties
        M;
        type = 'quad';
    end
    
    methods
        function obj = spec_quad(M, iwp, izp, rho)
            %SPEC_QUAD Construct an instance of this class
            %   Detailed explanation goes here
            if nargin > 3
                rho = 1;
            end            
            obj@spec_interface(iwp, izp, rho);
            obj.M = M;
        end
        
        function [M] = supply(obj)
            %SUPPLY quadratic performance specification
            M = obj.M;
        end
    end
end

