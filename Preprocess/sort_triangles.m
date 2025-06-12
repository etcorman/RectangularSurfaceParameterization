function [tri_ord,edge_ord,sign_edge] = sort_triangles(idx, T, E2T, T2T, E2V, T2E)

persistent tri_ord_cach edge_ord_cach sign_edge_cach;
if isempty(tri_ord_cach)
    nv = max(abs(T(:)));
    tri_ord_cach = cell(nv,1);
    edge_ord_cach = cell(nv,1);
    sign_edge_cach = cell(nv,1);
end

if isempty(tri_ord_cach{idx})
    [tri_ord,edge_ord,sign_edge] = sort_triangles_comp(idx, T, E2T, T2T, E2V, T2E);
    tri_ord_cach{idx} = tri_ord;
    edge_ord_cach{idx} = edge_ord;
    sign_edge_cach{idx} = sign_edge;
else
    tri_ord = tri_ord_cach{idx};
    edge_ord = edge_ord_cach{idx};
    sign_edge = sign_edge_cach{idx};
end