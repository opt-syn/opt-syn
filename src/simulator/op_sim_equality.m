classdef op_sim_equality < op_sim
    %OP_SIM_EQUALITY an operator used for the purposes of simulation (algorithm
    %execution). 

    %an affine operation enforcing an equality constraint Ez = b
    %
    %implemented as E'(Ez - b) = 0
    %with E full row rank

    %TODO: maybe this should be placed in the operator proper?
    
    properties        
        E= []; %forward evaluation (e.g. gradient)        
        b = [];  %function value (or function values in a game)        
        
    end
    
    methods
        function obj = op_sim_equality(E, b)
            %OP_SIM Construct an instance of this class
            %   operations used in the evaluation of the operator
            
            
            %create the operators
            %TODO: abstract to more general operators E (not stored in
            %memory as a matrix)

            fw = @(k, z, param) E' * (E*z - b);
            %https://proximity-operator.net/proximityoperator.html
            bw = @(k, D, z, param) z - E' * ((E*E') \ (E*z - b)); 
            f = @(k, z, param) norm(E*z - b);
            
            obj@op_sim(fw, bw, f);

            obj.E = E;
            if nargin < 2
                obj.b = 0;
            end
            obj.EQUALITY = 1;
        end               
    end
end

