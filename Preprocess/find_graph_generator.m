function [cycle,cocycle] = find_graph_generator(l, T, E2T, E2V, init)
% "Greedy Optimal Homotopy and Homology Generators"
% Jeff Erickson, Kim Whittlesey

if ~exist('init', 'var')
    init = 1;
end

nv = max(T(:));
ne = size(E2T,1);
nf = max(E2T(:));

%% Primal graph
% Set edge weights, zero at boundaries
ide_bound = any(E2T(:,1:2) == 0, 2);
l = abs(l) + 1e-5;
w = l;
w(ide_bound) = 0;

% Compute minimal spanning tree on primal graph
try
    G = graph(E2V(:,1), E2V(:,2), w);
    [Tree, pred] = minspantree(G, 'Type','tree', 'Root',init);
    pred = pred';
    EdgeTree = [pred(2:end), (2:nv)'];
catch
    warning('Kruskal algorithm from hell.');
    wAdj = sparse(E2V(:,[1 2]), E2V(:,[2 1]), [w,w], nv, nv);
    Adj = sparse(E2V(:,[1 2]), E2V(:,[2 1]), ones(nf,2), nv, nv);
    [~, EdgeTree, ~] = kruskal(Adj, wAdj);
    EdgeTree = sort(EdgeTree, 2);
    pred = tree_predecessor(init, EdgeTree);
end

% Find edge indices which are not in the tree
[~,id] = setdiff(sort(E2V,2), sort(EdgeTree,2), 'rows');
assert(all(~isnan(pred)), 'Multiple connected components.');

% Remove boundary edges
id = setdiff(id, find(ide_bound));

%% Dual graph
% Edge weights for the dual graph
w = 1./l;

% Compute minimal spanning tree on the dual graph of the remaining edges
try
    Gco = graph(E2T(id,1), E2T(id,2), w(id));
    [CoTree, copred] = minspantree(Gco, 'Type','forest', 'Method','sparse');
    copred = copred';
    EdgeCoTree = [copred(1:end), (1:max(vec(E2T(id,1:2))))'];%CoTree.Edges.EndNodes;
catch
    warning('Kruskal algorithm from hell.');
    wCoAdj = sparse(E2T(id,[1 2]), E2T(id,[2 1]), w, nf, nf);
    CoAdj = sparse(E2T(id,[1 2]), E2T(id,[2 1]), ones(length(id),2), nf, nf);
    [~, EdgeCoTree, ~] = kruskal(CoAdj, wCoAdj);
    copred = tree_predecessor(init, EdgeCoTree);
end

% Find edges indices which are not in the dual graph
[~,idCoEdge] = setdiff(sort(E2T(:,1:2),2), sort(EdgeCoTree,2), 'rows');

% Find edges which are neither in the primal nor in the dual graph
idGen = intersect(id, idCoEdge);

%% Build cycle basis
% Find the primal loops starting from idGen
cycle = cell(length(idGen),1);
for i = 1:length(idGen)
    left = flipud(predecessors(pred, E2V(idGen(i),1))); 
    right = predecessors(pred, E2V(idGen(i),2));
    cycle{i} = [left; right(1:end-1)];
end

% Find the dual loops starting from idGen
cocycle = cell(length(idGen),1);
assert(all(~isnan(pred)));
assert(all(~isnan(copred)));
for i = 1:length(idGen)
    left = flipud(predecessors(copred, E2T(idGen(i),1))); 
    right = predecessors(copred, E2T(idGen(i),2));
    cocycle{i} = [left; right(1:end-1)];
end

end

function path = predecessors(pred, i)
path = i;
i = pred(path(end));
while i ~= 0
    path = [path; i];
    i = pred(path(end));
end
end