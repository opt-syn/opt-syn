function M_sub = lmim_index(M, i, j)
%LMIM_INDEX subscript/index an LMIM matrix

[n, m] = dim(M);

if isempty(i)
    i = 1:n;
end
if isempty(j)
    i = 1:m;
end


li = length(i);
lj = length(j);
ei = sparse(1:li, i, ones(li, 1), li, n);
ej = sparse(1:lj, j, ones(lj, 1), lj, m)';


M_sub = full(ei) * M * full(ej);

end

