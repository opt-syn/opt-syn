classdef op_sim
    %OP_SIM an operator used for the purposes of simulation (algorithm
    %execution). 

    %TODO: maybe this should be placed in the operator proper?
    
    properties
        f;  %
        fw; %forward evaluation (e.g. gradient)
        bw; %backward evaluation (e.g. proximal operator)        
    end
    
    methods
        function obj = op_sim(f, fw, bw)
            %OP_SIM Construct an instance of this class
            %   operations used in the evaluation of the operator
            obj.f = f;
            obj.fw = fw;
            obj.bw = bw;
%             if nargin < 4
%                 obj.blocksize = [];
%             else
%                 obj.blocksize = blocksize;
%             end
        end               
    end
end

