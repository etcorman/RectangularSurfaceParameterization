function [idx_bound,edge_bound] = extract_all_boundary_curves(T, E2V)

deg = sum(T ~= 0, 2);
deg_max = max(deg);
deg_min = min(deg);
nf = size(T,1);

e2v = zeros(deg_max*nf,2);
for i = 1:deg_max
    for k = deg_min:deg_max
        j = mod(i, k) + 1;
        id = (i-1)*nf+1:i*nf;
        id = id(deg == k);
        e2v(id,:) = [T(deg == k,i), T(deg == k,j)];
    end
end
e2v = e2v(all(e2v ~= 0,2),:);
[~,ia,ic] = unique(sort(e2v,2),'rows');
E2V_u = e2v(ia,:);
n_adj = accumarray(ic, 1, [size(E2V_u,1),1]);
E2V_bound = E2V_u(n_adj == 1,:);

% e2v = [T(:,1), T(:,2); T(:,2), T(:,3); T(:,3), T(:,1)];
% [~,ia,ic] = unique(sort(e2v,2), 'rows');
% E2V_u = e2v(ia,:);
% n_adj = accumarray(ic, 1, [size(E2V_u,1),1]);
% E2V_bound = E2V_u(n_adj == 1,:);

if any(n_adj == 1)
    idx_bound = zeros(size(E2V_bound,1),1);
    id_start(1) = 1;
    k = 1;
    while ~isempty(E2V_bound)
        idx_bound(id_start(k):id_start(k)+1) = E2V_bound(1,:);
        E2V_bound(1,:) = [];

        for i = id_start(k)+2:size(idx_bound,1)+1
            id1 = E2V_bound(:,1) == idx_bound(i-1);
            id2 = E2V_bound(:,2) == idx_bound(i-1);
            if sum(id1) + sum(id2) <= 1
                if sum(id1) == 1
                    vx = E2V_bound(id1,2);
                    E2V_bound(id1,:) = [];
                elseif sum(id2) == 1
                    vx = E2V_bound(id2,1);
                    E2V_bound(id2,:) = [];
                else
                    error('Hole in boundary!?');
                end

                if vx ~= idx_bound(id_start(k))
                    idx_bound(i) = vx;
                else
                    id_start(k+1) = i;
                    k = k + 1;
                    break;
                end
            else
                error('Branching boundary!?'),
            end
        end
    end
%     id_start(end) = id_start(end)-1;
%     idx_bound = mat2cell(idx_bound, id_start(2:end),1);
    idx_bound = mat2cell(idx_bound, diff(id_start),1);
    
    n = length(idx_bound);
    edge_bound = cell(n,1);
    for i = 1:n
        [~,~,edge_bound{i}] = intersect(sort([idx_bound{i}, circshift(idx_bound{i}, [-1,0])], 2), E2V, 'rows', 'stable');
        assert(length(idx_bound{i}) == length(edge_bound{i}), 'Could not find all edges from boundary.');
    end
else
    idx_bound = [];
    edge_bound = [];
end
