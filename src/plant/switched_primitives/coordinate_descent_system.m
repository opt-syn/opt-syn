function [P] = coordinate_descent_system(c)
%COORDINATE_DESCENT_SYSTEM 
% block coordinate descent with c blocks
%
% can read out all subgradient information

Plist = coordinate_descent_primitives(c);
P = Plist{1};

end

