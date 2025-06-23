function [param,Src,dec] = preprocess_ortho_param(Src, dec, ifboundary, ifhardedge, tol_dihedral_deg, Ehard2V)
% Preprocess geometry:
% - Detect hard edges 
% - Compute boundary edges
% - Remesh so that each triangle has only one constrained edge
% - Store data for trivial connection (non-contractible cycles, Gaussian curvature, parallel transport)

% Input:
% - Src: triangle mesh data structure
% - dec: DEC data structure
% - ifboundary: (boolean) enforce boundary alignment
% - ifhardedge: (boolean) enforce hard-edge alignment
% - tol_dihedral_deg: (double) hard-edge detection threshold (angle in degree) 
% - Ehard2V: (integer array n x 2) vertex indices of alignment edges (optinal)

% Output:
% - param: data structure containing all parametrization constraints information
% - Src: remeshed triangle mesh data structure
% - dec: remeshed DEC data structure

comp_angle = @(u,v,n) atan2(dot(cross(u,v,2),n,2), dot(u,v,2));

%% Remeshing: a triangle cannot have two alignment constraints
% Check three cases:
% 1. alignment constraints are detected by the dihedral angle
% 2. alignment constraints is given by vertex indices

if ifhardedge && ~exist('Ehard2V','var')
    % Compute hard edges
    [ide_hard,tri_hard,ide_bound,tri_bound] = detect_hard_edge(Src, tol_dihedral_deg);
    
    % List of all constraints (boundary + hard edges)
    ide_fix = [ide_hard; ide_bound];    % Constrained edges
    tri_fix = [tri_hard(:); tri_bound]; % Corresponding faces
    
    % If a face appears twice it is over constrained
    % -> must be split in 3
    if numel(tri_fix) ~= length(unique(tri_fix))
        % Find indices of over-constrained triangles
        tri = sum(ismember(abs(Src.T2E), ide_fix) ,2) >= 2;

        % Add the barycenter to the new vertex list Xs
        b = (Src.X(Src.T(tri,1),:)+Src.X(Src.T(tri,2),:)+Src.X(Src.T(tri,3),:))/3;
        Xs = [Src.X; b];

        % Add 3 triangles per constrained triangles to the new triangle list Ts
        np = Src.nv+(1:size(b,1))';
        Ttri = [Src.T(tri,[1 2]), np ; Src.T(tri,[2 3]), np ; Src.T(tri,[3 1]), np];
        Ts = Src.T;
        Ts(tri,:) = [];
        Ts = [Ts; Ttri];

        % Recompute connectivity information
        Src = MeshInfo(Xs, Ts);

        % Recompute hard edges
        [ide_hard,tri_hard,ide_bound,tri_bound] = detect_hard_edge(Src, tol_dihedral_deg);
        
        % New list of all constraints (boundary + hard edges)
        ide_fix = [ide_hard; ide_bound];
        tri_fix = [tri_hard(:); tri_bound];
    end

    % Check that the remeshing worked
    assert(numel(tri_fix) == length(unique(tri_fix)), 'Multiple constraints on a triangle.');
elseif exist('Ehard2V','var') 
    [~,ide_hard] = intersect(Src.E2V, sort(Ehard2V,2), 'rows');
    tri_hard = Src.E2T(ide_hard,1:2);
else
    ide_hard = [];
    tri_hard = Src.E2T(ide_hard,1:2);
end

%% Remeshing: only constraint edges can have both vertices belonging to constraint edges or a boundary
% If this happen, the edge is plitted in two
ide_bound = boundary_indices(Src);
if ifhardedge || ~isempty(ide_bound)
    % Gather edge indices from constrained edges and boundary edges
    ide_fix = [ide_hard; ide_bound];
     % Find corresponding vertex indices
    idx_fix = unique(Src.E2V(ide_fix,:));

    % Edges indices whose vertices both belongs to constrained vertices
    ide = find(all(ismember(Src.E2V, idx_fix), 2));
    % Edges indices whose vertices both belongs to constrained vertices but
    % are not constrained edges themselves
    ide = setdiff(ide, ide_fix);

    % If ide is not empty, hese edges must be split
    if ~isempty(ide) 
        nv = Src.nv;            % Current vertex count
        Ts = Src.T;             % New face list
        Xs = Src.X;             % New vertex list
        E2T = Src.E2T(:,1:2);   % New edge to face table
        E2V = Src.E2V;          % New edge to vertices
        while ~isempty(ide) % while there are edges to split
            % Take first edge in the list
            id = ide(1);
            idx = E2V(id,:);
    
            % Split the first triangle in two
            idt1 = E2T(id,1);
            t1 = Ts(idt1,:);
            t1 = circshift(t1, [0,1-find(~ismember(t1,idx))]);
    
            % Split the second triangle in two
            idt2 = E2T(id,2);
            t2 = Ts(idt2,:);
            t2 = circshift(t2, [0,1-find(~ismember(t2,idx))]);
    
            % Update the triangle list
            Ts(idt1,:) = [t1(1), t1(2), nv+1];
            Ts(idt2,:) = [t2(1), t2(2), nv+1];
            Ts = [Ts; t1(1), nv+1, t1(3); t2(1), nv+1, t2(3)];
            
            % Middle of the edge is added to the vertex list
            m = (Xs(E2V(id,1),:) + Xs(E2V(id,2),:))/2;
            Xs = [Xs; m];
            nv = nv + 1;
    
            % Update connectivity and check for new edges to split
            idx = sort(E2V(ide_fix,:),2);
            [E2V,~,E2T] = connectivity(Ts);
            [~,ide_fix] = intersect(sort(E2V,2), idx, 'rows');
            idx_fix = unique(E2V(ide_fix,:));
            ide = find(all(ismember(E2V, idx_fix), 2));
            ide = setdiff(ide, ide_fix);
        end
    
        % Vertex pair corresponding to constrained edges (needed to find
        % new edge indices)
        idx_hard = sort(Src.E2V(ide_hard,:),2);

        % Recompute connectivity and DEC
        Src = MeshInfo(Xs, Ts);
        dec = dec_tri(Src);

        % Update constrained indices
        [~,ide_hard] = intersect(sort(Src.E2V,2), idx_hard, 'rows');
        tri_hard = Src.E2T(ide_hard,1:2);
    end
end

% Store hard edges
param.ide_hard = ide_hard;
param.tri_hard = tri_hard;

%% Compute boundary related stuff (mostly needed for, e.g., trivial connections)
[ide_bound,tri_bound] = boundary_indices(Src);
idx_bound = unique(Src.E2V(ide_bound,:));
idx_int = setdiff((1:Src.nv)', idx_bound);

ide_int = setdiff((1:Src.ne)', ide_bound);
tri_int = setdiff((1:Src.nf)', tri_bound);

% Store in structure
param.idx_bound = idx_bound;
param.ide_bound = ide_bound;
param.tri_bound = tri_bound;
param.idx_int = idx_int;
param.ide_int = ide_int;
param.tri_int = tri_int;

%% Merge boundary and hard edges
% group constraints by connected component  (still for trivial connections)
ide_fix = [ide_hard; ide_bound];

% Compute a graph made of only the constrained edges (still for trivial connections)
idx_fix = unique(Src.E2V(ide_fix,:));
idx_fix_inv = zeros(Src.nv,1);
idx_fix_inv(idx_fix) = 1:length(idx_fix);
G = graph(idx_fix_inv(Src.E2V(ide_fix,1)), idx_fix_inv(Src.E2V(ide_fix,2)));

% Find the connected components of the graph (still for trivial connections)
[bins,binsizes] = conncomp(G);

% For each component store: vertex, edge, triangle indices and edge
% orientation sign (still for trivial connections)
idx_fix_cell = cell(length(binsizes),1);
ide_fix_cell = cell(length(binsizes),1);
tri_fix_cell = cell(length(binsizes),1);
ide_sign_fix_cell = cell(length(binsizes),1);
for i = 1:length(binsizes)
    idx_fix_cell{i} = idx_fix(bins == i);
    ide_fix_cell{i} = find(all(ismember(Src.E2V, idx_fix_cell{i}), 2));
    ide_fix_cell{i} = ide_fix_cell{i}(ismember(ide_fix_cell{i}, ide_fix));
    tri_fix_cell{i} = vec(Src.E2T(ide_fix_cell{i},1:2));
    ide_sign_fix_cell{i} = vec(Src.E2T(ide_fix_cell{i},3:4));
end

%% Smooth cross field: Face-to-face parallel transport
E2T = zeros(Src.ne,2);
E2T(:,1) = Src.E2T(:,1).*(Src.E2T(:,3) > 0) + Src.E2T(:,2).*(Src.E2T(:,3) < 0);
E2T(:,2) = Src.E2T(:,1).*(Src.E2T(:,3) < 0) + Src.E2T(:,2).*(Src.E2T(:,3) > 0);

% Compute angle defect
K = gaussian_curvature(Src.X, Src.T);
assert(norm(sum(K) - 2*pi*(Src.nf-Src.ne+Src.nv)) < 1e-5, 'Gaussian curvature does not match topology.');

% Local basis: e1r aligned with constrained and boundary edges
edge = Src.X(Src.E2V(:,2),:) - Src.X(Src.E2V(:,1),:);
edge = edge./sqrt(sum(edge.^2,2));
e1r = Src.X(Src.T(:,2),:) - Src.X(Src.T(:,1),:);
e1r = e1r./repmat(sqrt(sum(e1r.^2, 2)), [1 3]);
for i = 1:length(binsizes)
    tri = tri_fix_cell{i};
    ide = [ide_fix_cell{i}; ide_fix_cell{i}];
    ides = ide_sign_fix_cell{i};
    e1r(tri(tri ~= 0),:) = edge(ide(tri ~= 0),:).*ides(tri ~= 0);
end
e2r = cross(Src.normal, e1r, 2);

% Angle between edge and local basis
edge_angles = zeros(Src.ne,2);
edge_angles(ide_int,1) = comp_angle(edge(ide_int,:), e1r(E2T(ide_int,1),:), Src.normal(E2T(ide_int,1),:));
edge_angles(ide_int,2) = comp_angle(edge(ide_int,:), e1r(E2T(ide_int,2),:), Src.normal(E2T(ide_int,2),:));

% Parallel transport
para_trans = wrapToPi(edge_angles(:,1) - edge_angles(:,2));
para_trans(ide_bound) = 0;

assert(norm(wrapToPi(dec.d1d*para_trans - K)) < 1e-6, 'Gaussian curvature incompatible with angle defect.');

% Angle between local basis and triangleedges
param.ang_basis = [comp_angle(Src.X(Src.T(:,1),:) - Src.X(Src.T(:,2),:), e1r, Src.normal), ...
                   comp_angle(Src.X(Src.T(:,2),:) - Src.X(Src.T(:,3),:), e1r, Src.normal), ...
                   comp_angle(Src.X(Src.T(:,3),:) - Src.X(Src.T(:,1),:), e1r, Src.normal)];

% Store stuff
param.E2T = E2T;
param.e1r = e1r;
param.e2r = e2r;
param.para_trans = para_trans;
param.Kt = K;
param.Kt_invisible = K - dec.d1d*para_trans;

%% Smooth cross field: take care of acute angle between edge constraints
% See: "Frame Fields for CAD models" https://inria.hal.science/hal-03537852/document

% Build exterior derivative of dual 1-form (d1d) where constrained edges
% are seen as boundaries
E2V = Src.E2V;
T = Src.T;
nv = Src.nv;
for i = idx_fix'
    [tri_ord,edge_ord,sign_edge] = sort_triangles(i, Src.T, Src.E2T, Src.T2T, Src.E2V, Src.T2E);
    id = ismember(edge_ord, ide_hard);
    n = sum(id);
    if (n == 1) && any(i == idx_int)
        continue;
    elseif (n == 1) && any(i == idx_bound)
        n = n + 1;
    end
    
    p = 1;
    flag = zeros(length(edge_ord),1);
    flag_tri = zeros(length(edge_ord),1);
    for j = 1:length(edge_ord)
        if id(j)
            p = mod(p, n) + 1;
        else
            flag(j) = p;
        end
        flag_tri(j) = p;
    end
    
    for j = 1:n-1
        ide = edge_ord(flag == j);
        E2V(ide,:) = (E2V(ide,:) ~= i).*E2V(ide,:) + (E2V(ide,:) == i).*(nv + j);

        tri = tri_ord(flag_tri == j);
        T(tri,:) = (T(tri,:) ~= i).*T(tri,:) + (T(tri,:) == i).*(nv + j);
    end

    if n > 1
        nv = nv + n - 1;
    end
end

ide_free = setdiff((1:Src.ne)', ide_fix);
d1d = sparse(E2V(ide_free,:), [ide_free, ide_free], [ones(length(ide_free),1),-ones(length(ide_free),1)], nv, Src.ne);

assert(all(sum(abs(d1d),2) ~= 0));

% Store stuff
Vp2V = unique([T(:), Src.T(:)], 'rows');
[~,id] = sort(Vp2V(:,1));
param.Vp2V = Vp2V(id,2);
param.d1d = d1d;
param.idx_fix_plus = [idx_fix; (Src.nv+1:nv)'];
param.idx_reg = setdiff((1:Src.nv)', idx_fix);

theta = angles_of_triangles(Src.X, Src.T);
K = 2*pi - accumarray(T(:), theta(:));
K(param.idx_fix_plus) = K(param.idx_fix_plus) - pi;
param.K = K;
param.K_invisible = K - param.d1d*para_trans;

%% Trivial connection: Path between isolated constraints
% Build the sparse matrix Ilink accumulating the dual 1-form along each
% path connecting the connected component of constraints
nc = max(length(ide_fix_cell) - 1, 0);
Ilink = sparse(nc,Src.ne);

% Compute a shortest path between the first component and all the other
for i = 1:nc
    % Edge weights
    ld = max(dec.star1p*sqrt(Src.SqEdgeLength), 1e-5);
    for j = 1:nc+1
        if (j ~= 1) && (j ~= i+1)
            ld(ide_fix_cell{j}) = max(ld)*1e5;  % Large weight on all boundaries except 1 and j
        else
            ld(ide_fix_cell{j}) = min(ld)*1e-5; % Small weight on all boundaries for 1 and j
        end
    end

    % Build prinmal graph
    Gd = graph(E2T(ide_int,1), E2T(ide_int,2), ld(ide_int));
    
    % Shortest dual path from a vertex in component 1 and a vertex in
    % compenent i+1
    P = shortestpath(Gd, tri_fix_cell{1}(1), tri_fix_cell{i+1}(1))';

    % Find edge indices of the path
    ed = [P(1:end-1), P(2:end)];
    [~,~,ide] = intersect(sort(ed,2), sort(E2T,2), 'rows', 'stable');
    assert(length(ide) == length(P)-1);
    
    % Remove part of the path on the boundary of 1 and i+1
    a = find(ismember(P, tri_fix_cell{1}), 1, 'last');
    b = find(ismember(P, tri_fix_cell{i+1}), 1, 'first');
    id = a:b-1;
    ide = ide(id);
    ed = ed(id,:);

    % Edge sign
    s = (E2T(ide,1) == ed(:,1)) - (E2T(ide,2) == ed(:,1));

    % Build 1-form integration matrix along the path
    Ilink(i,:) = sparse(ones(length(ide),1), ide, s, 1, Src.ne);
end
Ilink_hard = Ilink;
Ilink_hard(:,ide_hard) = 0; % hard edge do not count

param.Ilink = Ilink;
param.Ilink_hard = Ilink_hard;

%% Trivial connection: Non-contractible cycles
% Compute non-constractible loops
[cycle,cocycle] = find_graph_generator(full(diag(dec.star1p)), Src.T, Src.E2T, Src.E2V, 1);

% Build the sparse matrix Icycle accumulating the dual 1-form along each
% cycle
nc = length(cocycle);
Icycle = sparse(nc,Src.ne);
for i = 1:nc
    % Find edge indices correpsonding to cycle
    ed = [cocycle{i}, circshift(cocycle{i}, [1,0])];
    [~,ide] = ismember(sort(ed,2), sort(E2T,2), 'rows');
    assert(length(ide) == length(cocycle{i}));
    
    % Find edge sign in the cycle
    s = (E2T(ide,1) == ed(:,1)) - (E2T(ide,2) == ed(:,1));

    % Build sparce matrix
    Icycle(i,:) = sparse(ones(length(ide),1), ide, s, 1, Src.ne);
end
Icycle_hard = Icycle;
Icycle_hard(:,ide_hard) = 0; % hard edge do not count

param.Icycle = Icycle;
param.Icycle_hard = Icycle_hard;

%% Set constraint list
if ifboundary && ifhardedge
    ide_fix = [ide_hard; ide_bound];
    tri_fix = [tri_hard(:); tri_bound];
elseif ifboundary
    ide_fix = ide_bound;
    tri_fix = tri_bound;
elseif ifhardedge
    ide_fix = ide_hard;
    tri_fix = tri_hard(:);
else
    ide_fix = [];
    tri_fix = [];
end

idx_fix = unique(Src.E2V(ide_fix,:));
param.ide_fix = ide_fix;
param.idx_fix = idx_fix;
param.ide_free = setdiff((1:Src.ne)', ide_fix);
param.tri_fix = tri_fix;
param.tri_free = setdiff((1:Src.nf)', param.tri_fix);
end

function [ide_hard,tri_hard,ide_bound,tri_bound] = detect_hard_edge(Src, tol_dihedral_deg)
comp_angle = @(u,v,n) atan2(dot(cross(u,v,2),n,2), dot(u,v,2));
tol_dihedral = tol_dihedral_deg*pi/180;

% Find boundary edges
[ide_bound,tri_bound] = boundary_indices(Src);

% Interior edges
ide_int = setdiff((1:Src.ne)', ide_bound);

% Compute unit vector edges
edge = Src.X(Src.E2V(:,2),:) - Src.X(Src.E2V(:,1),:);
edge = edge./sqrt(sum(edge.^2,2));

% Compute angle between normals by using the edge vector for sign
dihedral_angle = Src.E2T(ide_int,4).*comp_angle(Src.normal(Src.E2T(ide_int,1),:), Src.normal(Src.E2T(ide_int,2),:), edge(ide_int,:));

% Hard edges are angle larger than a threshold
ide_hard = ide_int(abs(dihedral_angle) > tol_dihedral);

% Find adjancent triangles
tri_hard = Src.E2T(ide_hard,1:2);
end

function [ide_bound,tri_bound] = boundary_indices(Src)
% A boundary edge belongs to only one triangle
ide_bound = find(any(Src.E2T == 0, 2));

% A boundary triangle is incident to a boundary edge
tri_bound = sum(Src.E2T(ide_bound,1:2),2);
end
