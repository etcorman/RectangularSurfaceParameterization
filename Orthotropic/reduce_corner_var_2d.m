function [Edge_jump,v2t,base_tri] = reduce_corner_var_2d(Src)

base_tri = zeros(Src.nv,1);
path_vx = cell(Src.nv,1);
path_edge = cell(Src.nv,1);
Tc = reshape((1:3*Src.nf)', [Src.nf,3]);
for i = 1:Src.nv
    [tri_ord,edge_ord,sign_edge] = sort_triangles(i, Src.T, Src.E2T, Src.T2T, Src.E2V, Src.T2E);
    edge_ord = edge_ord.*sign_edge;
    
    base_tri(i) = tri_ord(1);
    
    if tri_ord(end) == 0
        idvx = sum(Tc(tri_ord(1:end-1),:) .* (Src.T(tri_ord(1:end-1),:) == i), 2);
        path_vx{i} = [idvx, i*ones(length(idvx),1)];
        
        I = repelem(idvx(2:end), (1:length(idvx)-1)');
        J = zeros(size(I));
        d = edge_ord(2:end-1);
        id = 1;
        for k = 1:length(idvx)-1
            J(id) = d(1:k);
            id = id(end) + (1:k+1);
        end
        path_edge{i} = [I, J];
    else
        idvx = sum(Tc(tri_ord,:) .* (Src.T(tri_ord,:) == i), 2);
        path_vx{i} = [idvx, i*ones(length(idvx),1)];
        
        I = repelem(idvx(2:end), (1:length(idvx)-1)');
        J = zeros(size(I));
        d = edge_ord(2:end);
        id = 1;
        for k = 1:length(idvx)-1
            J(id) = d(1:k);
            id = id(end) + (1:k+1);
        end
        path_edge{i} = [I, J];
    end
end

I = cell2mat(path_edge);
Edge_jump = sparse(I(:,1), abs(I(:,2)), sign(I(:,2)), 3*Src.nf, Src.ne);

I = cell2mat(path_vx);
v2t = sparse(I(:,1), I(:,2), 1, 3*Src.nf, Src.nv);
