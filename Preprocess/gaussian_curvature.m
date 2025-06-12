function [K,idx_bound,ide_bound] = gaussian_curvature(X, T)
% Output gaussian curvature inside and normal curvature at the boundary

% Compute interior curvature
theta = angles_of_triangles(X, T);
K = 2*pi - accumarray(T(:), theta(:));
% K(idx_bound) = K(idx_bound) - pi;

% Compute boundary normal curvature
E2V = sort([T(:,1), T(:,2); T(:,2), T(:,3); T(:,3), T(:,1)], 2);
[E2V_u,~,ic] = unique(E2V, 'rows');
n_adj = accumarray(ic, 1, [size(E2V_u,1),1]);
ide_bound = find(n_adj == 1);
E2V_bound = E2V_u(ide_bound,:);
idx_bound = unique(E2V_bound);
K(idx_bound) = K(idx_bound) - pi;
