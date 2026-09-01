function [s_array] = lmim_vec_break(s)
%LMIM_VEC_BREAK break up an lmim array into a vector of lmim objects
sz = ssize(s, 1);

% I = speye(sz);
s_array = lmim.empty(sz, 0);
for i = 1:sz
    e = zeros(1, sz);
    e(i) = 1;
    s_array(i) = e*s;
end

end

