function [M_tilde] = loop_flip(M)
%LOOP_FLIP define the LFT factor for the static filter
%Convert [p; q] = M [z w] to [p; w] = Mtilde [z; q]
%used for passivity conversion

assert(M(2, 2) ~= 0)

drec = 1/M(2, 2);
M_tilde = [M(1, 1) - M(1, 2)*M(2, 1)*drec, M(1, 2)*drec;
    -M(2, 1)*drec, drec];



end

