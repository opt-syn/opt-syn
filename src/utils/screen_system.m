function [E_indexer] = screen_system(n, ind)
%SCREEN_SYSTEM index the channels of a partitioned linear system (ss)
%  
%Args
%   n (array of int):              dimensions of partition
%   ind (cell of array of int):    indices 


N = sum(n);
k = length(n);

l = cellfun(@length, ind);

nc = cumsum([0, n]);
n_ind = sum(l);

i = zeros(1, n_ind);
j = zeros(1, n_ind);
v = ones(1, n_ind);


count = 0;
for c = 1:k
    na_curr = length(ind{c});

    i(count + (1:na_curr)) = count + (1:na_curr);
    j(count + (1:na_curr)) = nc(c) + ind{c};   

    count = count + na_curr;
end


E_indexer = full(sparse(i, j, v, n_ind, N));


end