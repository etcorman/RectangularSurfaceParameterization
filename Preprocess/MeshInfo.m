function mesh = MeshInfo(X, T)

mesh.X = X;
mesh.T = T;
assert(size(T,2) == 3, 'Not a triangulations.');
[mesh.E2V, mesh.T2E, mesh.E2T, mesh.T2T] = connectivity(mesh.T);

mesh.nf = size(mesh.T,1);
mesh.nv = size(mesh.X,1);
mesh.ne = size(mesh.E2V,1);

% Normals and areas
mesh.normal = cross(mesh.X(mesh.T(:,1),:) - mesh.X(mesh.T(:,2),:), mesh.X(mesh.T(:,1),:) - mesh.X(mesh.T(:,3),:));
mesh.area = sqrt(sum(mesh.normal.^2, 2))/2;
mesh.normal = mesh.normal./repmat(sqrt(sum(mesh.normal.^2, 2)), [1, 3]);

A = sparse(mesh.T, repmat((1:mesh.nf)', [3,1]), repmat(mesh.area, [3,1]), mesh.nv, mesh.nf);
mesh.Nv = A*mesh.normal;
mesh.Nv = mesh.Nv./repmat(sqrt(sum(mesh.Nv.^2,2)), [1,3]);

% Edge length
mesh.SqEdgeLength = sum((mesh.X(mesh.E2V(:,1),:) - mesh.X(mesh.E2V(:,2),:)).^2, 2);

% Angles
mesh.corner_angle = angles_of_triangles(mesh.X, mesh.T);
mesh.cot_corner_angle = cot(mesh.corner_angle);

end
