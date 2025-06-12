function [tri_ord,edge_ord,sign_edge] = sort_triangles_comp(idx, T, E2T, T2T, E2V, T2E)
% sort ring triangles around idx
ifbound = any(vec(E2V(any(E2T(:,1:2) == 0, 2),:) == idx));

if ifbound
    edi = E2T(any(E2V == idx, 2),:);
    
    m = edi(:,1:2) == 0;
    
    vxi = E2V(any(E2V == idx, 2),:);
    vxi = vxi(any(m ~= 0,2),:);
    vxi = sum(vxi .* (vxi ~= idx),2);
    vxi = [idx*ones(2,1), vxi];
    
    tri_start = edi(:,1:2).*edi(:,3:4).*m(:,[2 1]);
    tri_start = abs(tri_start(tri_start ~= 0));
    
    if1 = false; deg1 = sum(T(tri_start(1),:) ~= 0);
    if2 = false; deg2 = sum(T(tri_start(2),:) ~= 0);
    for k = 1:size(T,2)
        i = mod(k-1, deg1) + 1;
        j = mod(k  , deg1) + 1;
        if1 = if1 || all(T(tri_start(1),[i j]) == vxi(1,:));
        i = mod(k-1, deg2) + 1;
        j = mod(k  , deg2) + 1;
        if2 = if2 || all(T(tri_start(2),[i j]) == vxi(2,:));
    end
    
    if if1
        tri_start = tri_start(1);
    elseif if2
        tri_start = tri_start(2);
    else
        error('Could not find a staring triangle.');
    end
    
%     if all(T(tri_start(1),[1 2]) == vxi(1,:)) || all(T(tri_start(1),[2 3]) == vxi(1,:)) || all(T(tri_start(1),[3 1]) == vxi(1,:))
%         tri_start = tri_start(1);
%     elseif all(T(tri_start(2),[1 2]) == vxi(2,:)) || all(T(tri_start(2),[2 3]) == vxi(2,:)) || all(T(tri_start(2),[3 1]) == vxi(2,:))
%         tri_start = tri_start(2);
%     else
%         error('Could not find a staring triangle.');
%     end
else
    ifbound = false;
end

trii = find(any(T == idx, 2));
if ifbound
    trii = circshift(trii, [1-find(trii == tri_start),0]);
end
tri_ord = zeros(length(trii),1);
tri_ord(1) = trii(1);
tovisit = trii(2:end);
k = 1;
while ~isempty(tovisit)
    [t,~,id] = intersect(T2T(tri_ord(k),:), tovisit);
    
    % pick a consistent ordering
    if length(t) == 2
        if k == 1
            t1 = T(t(1),:);
            t1 = circshift(t1, [0,1-find(t1 == idx)]);
            
            if all(ismember(t1(1:2), T(trii(1),:)))
                t = t(1);
                id = id(1);
            else
                t = t(2);
                id = id(2);
            end
        else
            error('Problem when sorting triangles.');
        end
    end
    
    k = k + 1;
    tri_ord(k) = t;
    tovisit(id) = [];
end
if ifbound
    tri_ord = [tri_ord; 0];
end

% sort edges
if nargout >= 2 
    if length(tri_ord) > 2
        path = [tri_ord, circshift(tri_ord, [1,0])];
        paths = sort(path,2);
        if exist('T2E','var')
            ide = unique(abs(T2E(trii,:)));
            [~,~,edge_ord] = intersect(paths, sort(E2T(ide,1:2),2), 'rows', 'stable');
            edge_ord = ide(edge_ord);
        else
            [~,~,edge_ord] = intersect(paths, sort(E2T(:,1:2),2), 'rows', 'stable');
        end
        sign_edge = (path(:,1) == E2T(edge_ord,1)) .* E2T(edge_ord,3) + (path(:,1) == E2T(edge_ord,2)) .* E2T(edge_ord,4);
        
        if length(edge_ord) ~= length(tri_ord)
            assert(size(unique(paths, 'rows'),1) == size(paths,1), 'Non manifold mesh.');
            error('Could not find edge ordering.');
        end
    else
        path = [tri_ord, circshift(tri_ord, [1,0])];
        edge_ord = find((E2T(:,1) == tri_ord(1)) & (E2T(:,2) == tri_ord(2)));
        sign_edge = (path(:,1) == E2T(edge_ord,1)) .* E2T(edge_ord,3) + (path(:,1) == E2T(edge_ord,2)) .* E2T(edge_ord,4);
    end
    
end