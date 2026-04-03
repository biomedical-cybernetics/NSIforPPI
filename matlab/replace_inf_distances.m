function replaced_matrix = replace_inf_distances(distance_matrix)
inf_indices = isinf(distance_matrix);
if ~any(inf_indices(:))
    replaced_matrix = distance_matrix;
    return;
end
G = graph(distance_matrix ~= Inf);
components = conncomp(G);
max_distances = zeros(1, max(components));
for i = 1:max(components)
    component_indices = components == i;
    component_distances = distance_matrix(component_indices, component_indices);
    max_distances(i) = max(component_distances(:));
end
replaced_matrix = distance_matrix;
replaced_matrix(inf_indices) = sum(max_distances);
end
