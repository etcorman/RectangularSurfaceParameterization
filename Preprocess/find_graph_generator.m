function [cycle,cocycle] = find_graph_generator(X, l, T, E2T, E2V, init)

if ~exist('init', 'var')
    init = 1;
end

nv = max(T(:));
ne = size(E2T,1);
nf = max(E2T(:));

ide_bound = any(E2T(:,1:2) == 0, 2);

l = abs(l) + 1e-5;
% w = [sqrt(sum((X(T(:,1),:)-X(T(:,2),:)).^2,2)), sqrt(sum((X(T(:,2),:)-X(T(:,3),:)).^2,2)), sqrt(sum((X(T(:,3),:)-X(T(:,1),:)).^2,2))];
w = l;
w(ide_bound) = 0;
try
    G = graph(E2V(:,1), E2V(:,2), w);
    [Tree, pred] = minspantree(G, 'Type','tree', 'Root',init);
    pred = pred';
    EdgeTree = [pred(2:end), (2:nv)'];%Tree.Edges.EndNodes;
catch
    warning('Kruskal algorithm from hell.');
    wAdj = sparse(E2V(:,[1 2]), E2V(:,[2 1]), [w,w], nv, nv);
    Adj = sparse(E2V(:,[1 2]), E2V(:,[2 1]), ones(nf,2), nv, nv);
%     wAdj = sparse(T(:,[1 2 3]), T(:,[2 3 1]), w, nv, nv);
%     Adj = sparse(T(:,[1 2 3]), T(:,[2 3 1]), ones(nf,3), nv, nv);
    [~, EdgeTree, ~] = kruskal(Adj, wAdj);
    EdgeTree = sort(EdgeTree, 2);
    pred = tree_predecessor(init, EdgeTree);
end
[~,id] = setdiff(sort(E2V,2), sort(EdgeTree,2), 'rows');
assert(all(~isnan(pred)), 'Multiple connected components.');

id = setdiff(id, find(ide_bound));
bar = (X(T(:,1),:) + X(T(:,2),:) + X(T(:,3),:))/3;
mid = (X(E2V(:,1),:) + X(E2V(:,2),:))/2;
ld = zeros(ne,1);
ld(E2T(:,1) ~= 0) = sqrt(sum((mid(E2T(:,1) ~= 0) - bar(E2T(E2T(:,1) ~= 0,1))).^2,2));
ld(E2T(:,2) ~= 0) = ld(E2T(:,2) ~= 0) + sqrt(sum((mid(E2T(:,2) ~= 0) - bar(E2T(E2T(:,2) ~= 0,2))).^2,2));
w = -[ld(id),ld(id)]; w = w -min(w(:))+0.1;
w = 1./l;
assert(all(w(:,1) > 0))
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
[~,idCoEdge] = setdiff(sort(E2T(:,1:2),2), sort(EdgeCoTree,2), 'rows');

idGen = intersect(id, idCoEdge);
cycle = cell(length(idGen),1);
for i = 1:length(idGen)
    left = flipud(predecessors(pred, E2V(idGen(i),1))); 
    right = predecessors(pred, E2V(idGen(i),2));
    cycle{i} = [left; right(1:end-1)];
end

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