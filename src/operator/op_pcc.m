classdef op_pcc < op_sml
    %OP_PCC An operator which is the subdifferential of a proper, closed,
    %convex function.
    
    % Example: F = partial I_K, where I_K is the indicator function of a
    % closed convex set K.
    %
    %
    % noncausal multipliers    



    methods
        function obj = op_pcc(c)
            %OP_PCC Constructor for op_sml(0, inf, c)
            %
            %Args:
            %   c:  coordinate dimension

            if nargin < 1
                c = 1;
            end
            obj@op_sml(0 ,inf , c)      

        end

    end
end

