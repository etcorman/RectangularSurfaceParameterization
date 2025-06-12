function [A] = angles_of_triangles(V, T)
    % Computes for each triangle the 3 angles among its edges.
    % Input:
    %   option 1:   [A] = angles_of_triangles(V, T)
    %
    %               V   - (num_of_vertices x 3) 3D coordinates of
    %                     the mesh vertices.
    %               T   - (num_of_triangles x 3) T[i] are the 3 indices
    %                     corresponding to the 3 vertices of the i-th
    %                     triangle. The indexing is based on -V-.
    %
    %   option 2:   [A] = angles_of_triangles(L)
    %
    %               L - (num_of_triangles x 3) L[i] is a triple
    %                   containing the lengths of the 3 edges
    %                   corresponding to the i-th triange.
    %
    % Output:
    %

	E1 = V(T(:,2),:) - V(T(:,1),:);
    E2 = V(T(:,3),:) - V(T(:,2),:);
    E3 = V(T(:,1),:) - V(T(:,3),:);

    E1 = E1./sqrt(sum(E1.^2,2));
    E2 = E2./sqrt(sum(E2.^2,2));
    E3 = E3./sqrt(sum(E3.^2,2));
    A = pi - acos([dot(E3, E1, 2), dot(E1, E2, 2), dot(E2, E3, 2)]);
    assert(all(~isnan(A(:))) && all(isreal(A(:))), 'Triangle of size zero.');
end