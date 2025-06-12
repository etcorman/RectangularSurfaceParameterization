function [Edge_jump,v2t,base_tri] = reduce_corner_var_2d_cut(Src, ide_cut)

if ~exist('ide_cut','var') || isempty(ide_cut)
    ide_cut = false(Src.ne,1);
end
if ~islogical(ide_cut)
    ide_cutb = false(Src.ne,1);
    ide_cutb(ide_cut) = true;
    ide_cut = ide_cutb;
end

% Reduce variables at corner
base_tri = zeros(Src.nf,3);
path_vx = cell(Src.nv,1);
path_edge = cell(Src.nv,1);
Tc = reshape((1:3*Src.nf)', [Src.nf,3]);
nv = Src.nv;
for i = 1:Src.nv
    [tri_ord,edge_ord,sign_edge] = sort_triangles(i, Src.T, Src.E2T, Src.T2T, Src.E2V, Src.T2E);
    edge_ord = edge_ord.*sign_edge;
    if tri_ord(end) == 0
        ifbound = true;
        
        ide_cut(abs(edge_ord([1 end]))) = false;
        tri_ord(end) = [];
        edge_ord(end) = [];
    else
        ifbound = false;
        
        id = ide_cut(abs(edge_ord));
        if any(id)
            tri_ord = circshift(tri_ord, [1-find(id, 1, 'first'),0]);
            edge_ord = circshift(edge_ord, [1-find(id, 1, 'first'),0]);
            
            id = ide_cut(abs(edge_ord));
            assert(id(1))
        end
    end
    id = ide_cut(abs(edge_ord));
    n = sum(id);
    
    p = 0;
    flag = zeros(length(edge_ord),1);
    for j = 1:length(edge_ord)
        if id(j)
            p = p + 1;
        end
        flag(j) = p;
    end
    if ~ifbound
        flag(flag == n) = 0;
        n = max(0,n-1);
    end
    
    for j = 0:n
        % equivalence between vertices
        idt = tri_ord(flag == j);
        idvx = sum(Tc(idt,:) .* (Src.T(idt,:) == i), 2);
        assert(~isempty(idvx))
        if j == 0
            path_vx{i} = [idvx, i*ones(length(idvx),1)];
        else
            path_vx{nv+j} = [idvx, (nv+j)*ones(length(idvx),1)];
        end
        
        % edge path linking vertices
        I = repelem(idvx(2:end), (1:length(idvx)-1)');
        J = zeros(size(I));
        d = edge_ord(flag == j);
        d(1) = [];
        id = 1;
        for k = 1:length(idvx)-1
            J(id) = d(1:k);
            id = id(end) + (1:k+1);
        end
        if j == 0
            path_edge{i} = [I, J];
            if isempty(path_edge{i})
                path_edge{i} = zeros(0,2);
            end
        else
            path_edge{nv+j} = [I, J];
            if isempty(path_edge{nv+j})
                path_edge{nv+j} = zeros(0,2);
            end
        end
        
        % base triangle
        base_tri(idvx) = idt(1);
    end
    if n >= 1
        nv = nv + n;
    end
end

I = cell2mat(path_edge);
Edge_jump = sparse(I(:,1), abs(I(:,2)), sign(I(:,2)), 3*Src.nf, Src.ne);

I = cell2mat(path_vx);
v2t = sparse(I(:,1), I(:,2), 1, 3*Src.nf, nv);
