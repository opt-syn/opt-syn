function walks = all_walks(G, n, curr_node, curr_walk)
%ALL_WALKS: find all walks of length n in the graph G starting from the node curr_node

	if n == 0
		walks = curr_walk;
		return;
	end

	walks = [];
	nb = find(G(curr_node, :));

	for i = 1:length(nb)
		next_node = nb(i);
		next_walk = [curr_walk, next_node];
		walk_sub = all_walks(G, n-1, next_node, next_walk);
		walks = [walks; walk_sub];
	end	

end