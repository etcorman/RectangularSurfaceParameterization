function [SrcCut,idx_cut_inv,ide_cut_inv,edge_cut] = cut_mesh(X, T, E2V, E2T, T2E, T2T, idcone, edge_cut_tag)

nv = size(X,1);
ne = size(E2V,1);
nf = size(T,1);


% compute tree from first cone or boundary
[idx_bound_cell,ide_bound_cell] = extract_all_boundary_curves(T, E2V);
idx_bound = cell2mat(idx_bound_cell);
ide_bound = cell2mat(ide_bound_cell);

% dual spanning tree
Q = 1;
pred =-ones(nf,1); % Visited vertices
pred(Q) = 0;
edge = false(ne,1);
edge(edge_cut_tag) = true;
while ~isempty(Q)
    idtri = Q(1);
    Q(1) = [];

    adj = T2T(idtri,:);
    adj(adj == 0) = [];
    adjedge = abs(T2E(idtri,:));
    [~,ib] = ismember(adj, sum(E2T(adjedge,1:2).*(E2T(adjedge,1:2) ~= idtri),2));
    adjedge = adjedge(ib);
    
    for i = 1:length(adj)
        if (pred(adj(i)) == -1) && ~edge(adjedge(i))
            pred(adj(i)) = idtri;
            Q = [Q; adj(i)];
            
            edge(adjedge(i)) = true;
        end
    end
end

% cut leaves (ie degree vertices that are not singularities)
edge(edge_cut_tag) = false; 
edge_cut =~edge; % Set of current edge in the cut
deg = accumarray(vec(E2V(edge_cut,:)), 1, [nv,1]); % vertices degree

deg1 = deg == 1; % set of degree 1 vertices
while sum(deg1(idcone)) ~= sum(deg1) % while there's degree 1 vertices
    deg1(idcone(deg1(idcone))) = false;
    edge_cut(any(deg1(E2V),2)) = false; % remove edges attached to degree 1 vertices from the cut
    deg = accumarray(vec(E2V(edge_cut,:)), 1, [nv,1]); % recompute vertices degree with the new cut
    
    deg1 = deg == 1; % set of degree 1 vertices
end

% Create new mesh with cut edge as boundary
if any(edge_cut) && any(all(E2T(edge_cut,1:2) ~= 0,2))
    % Compute equivalence based on edges
    e2v = [T(:,1),T(:,2); T(:,2),T(:,3); T(:,3),T(:,1)];
    [e2vs,ids] = sort(e2v,2);
    id_cut = find(ismember(e2vs, E2V(edge_cut,:), 'rows'));
    
    [~,ia,ic] = unique(e2vs, 'rows');
    idtri_bound = ia(accumarray(ic,1) == 1);
    idtri1 = setdiff(ia(accumarray(ic,1) == 2), id_cut);
    idtri2 = setdiff((1:size(e2vs,1))', [idtri1; idtri_bound; id_cut]);
    [~,id1,id2] = intersect(e2vs(idtri1,:), e2vs(idtri2,:), 'rows');
    
    Tc = reshape((1:3*nf)', [nf,3]);
    e2v_cut = [Tc(:,1),Tc(:,2); Tc(:,2),Tc(:,3); Tc(:,3),Tc(:,1)];
    e2v_cuts = [e2v_cut(:,1).*(ids(:,1) == 1) + e2v_cut(:,2).*(ids(:,1) == 2), e2v_cut(:,1).*(ids(:,2) == 1) + e2v_cut(:,2).*(ids(:,2) == 2)];
    equiv_vx = [e2v_cuts(idtri1(id1),1), e2v_cuts(idtri2(id2),1) ; e2v_cuts(idtri1(id1),2), e2v_cuts(idtri2(id2),2)];
    
    % Propagate equivalences
    Tc = reshape(union_find(3*nf, equiv_vx), [nf,3]);
    idx_cut_inv = unique([Tc(:), T(:)], 'rows');
    idx_cut_inv = idx_cut_inv(:,2);
    assert(max(Tc(:)) == length(idx_cut_inv), 'Failure to find new indices.');
    
    Xc = X(idx_cut_inv,:);
else
    idx_cut_inv = (1:nv)';
    Tc = T;
    Xc = X;
end

SrcCut = MeshInfo(Xc, Tc);
chiCut = SrcCut.nf - SrcCut.ne + SrcCut.nv;
if chiCut ~= 1
    warning('Not topological disk after cut.');
end

[~,ide_cut_inv] = ismember(sort(idx_cut_inv(SrcCut.E2V),2), E2V, 'rows');
assert(all(ide_cut_inv ~= 0))
ids = idx_cut_inv(SrcCut.E2V(:,1)) == E2V(ide_cut_inv,1);
ide_cut_inv = ids.*ide_cut_inv - (~ids).*ide_cut_inv;

% figure;
% trisurf(SrcCut.T, SrcCut.X(:,1), SrcCut.X(:,2), SrcCut.X(:,3));
% axis equal;
end

function x = union_find(n, equiv)
assert(size(equiv,2) == 1 || size(equiv,2) == 2, 'Argument size invalid.');
assert(all(equiv(:) > 0), 'Wrong indexes.');

parent = (1:n)';
if size(equiv,2) == 2
    assert(all(equiv(:) <= n), 'Wrong indexes.');
    
    for i = 1:size(equiv,1)
        parent = union_tree(equiv(i,1), equiv(i,2), parent);
    end
elseif size(equiv,2) == 1
    neq = max(equiv);
    for i = 1:neq
        ind = mod(find(equiv == i)-1,n)+1;
        nind = length(ind);
        for j = 1:nind
            parent = union_tree(ind(j), ind(mod(j,nind)+1), parent);
        end
    end
else
    error('Argument size invalid.');
end

% Set nodes to root
x = parent;
for i = 1:n
    [x(i),~] = find_root(i, parent);
end

% Rearrange root nodes
nset = unique(x);
uniq_id = zeros(max(nset),1);
uniq_id(nset) = 1:length(nset);
x = uniq_id(x);
end

function parent = union_tree(x, y, parent)
[x_root,x_size] = find_root(x, parent);
[y_root,y_size] = find_root(y, parent);

if x_root ~= y_root
    if x_size > y_size
        parent(y_root) = x_root;
    else
        parent(x_root) = y_root;
    end
end
end

function [root,s] = find_root(x, parent)
root = x;
s = 0;
while root ~= parent(root)
    root = parent(root);
    s = s + 1;
end
end
