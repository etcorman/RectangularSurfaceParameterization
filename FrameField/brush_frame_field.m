function [ang,init_tree] = brush_frame_field(param, omega, tri_fix, ang_init)

if ~exist('ang_init', 'var')
    ang_init = zeros(length(tri_fix),1);
end

if ~isempty(tri_fix)
    init_tree = tri_fix(1);
elseif ~isempty(param.tri_bound)
    init_tree = param.tri_bound(1);
else
    init_tree = 1;
end

nf = max(param.E2T(:));
ang = zeros(nf,1);
ang(tri_fix) = ang_init;
ang = breadth_first_search(ang, omega(param.ide_int) - param.para_trans(param.ide_int), param.E2T(param.ide_int,:), @(x,y) x+y, init_tree);

end

function [y,S,level] = breadth_first_search(x, omega, E2V, fun, init)
% Input:
% x     : variable defined at vertices
% omega : define edges
% E2V   : graph edges
% fun   : update function
%
% Output:
% y: defined vertices, propagation of x(1) in the graph updated at edges
% with the function: x_i = fun(x_j,omega_ij)
% S: spanning tree creating during bread first search

if ~exist('init', 'var')
    init = 1;
end

nv = max(E2V(:));
ne = size(E2V,1);

assert(size(x,1) == nv, 'Variable has wrong dimension.');
assert(size(omega,1) == ne, 'Update has wrong dimension.');
y = x;

Q = init;
S = -ones(nv,1); % Visited vertices
S(Q) = 0;
level = zeros(nv,1);
l = 0;
while ~isempty(Q)
    idx = Q(1);
    Q(1) = [];
    l = l + 1;
    
    id1 = find(E2V(:,1) == idx);
    id2 = find(E2V(:,2) == idx);
    adjedge = [id1,-ones(size(id1)); id2, ones(size(id2))];
    adj = [E2V(id1,2)', E2V(id2,1)'];
    adjedge(adj == 0,:) = [];
    adj(adj == 0) = [];
    for i = 1:length(adj)
        if S(adj(i)) == -1
            S(adj(i)) = idx;
            level(adj(i)) = l;
            Q = [Q; adj(i)];
            
            s = adjedge(i,2);
            y(adj(i),:) = fun(y(idx,:), s*omega(adjedge(i,1),:));
        end
    end      
end
end