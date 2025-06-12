function [SrcCut,dec_cut,Align,Rot] = mesh_to_disk_seamless(Src, param, ang, sing, k21, ifseamless_const, ifboundary, ifhardedge)

% Cut geometry
idcone = param.idx_int(abs(sing(param.idx_int)) > 0.1); % Cone indices
[SrcCut,~,ide_cut_inv,~] = cut_mesh(Src.X, Src.T, Src.E2V, Src.E2T, Src.T2E, Src.T2T, idcone, k21 ~= 1);
dec_cut = dec_tri(SrcCut);

% Build hard seamless constraints
if ifseamless_const && nargout > 2
    % Find equivalences between edges in the cut and the boundaries of 'SrcCut'
    edge_bound_cut = find(any(SrcCut.E2T(:,1:2) == 0, 2));
    edge_bound_cut = [edge_bound_cut, sum(SrcCut.E2T(edge_bound_cut,1:2),2)];
    edge_bound_cut = edge_bound_cut(~ismember(abs(ide_cut_inv(edge_bound_cut(:,1))), param.ide_bound),:);
    [~,id] = sort(abs(ide_cut_inv(edge_bound_cut(:,1))));
    ide_cut = ide_cut_inv(edge_bound_cut(id,1));
    assert(all(abs(ide_cut(1:2:end-1)) == abs(ide_cut(2:2:end))), 'Cut failed.');
    ide_cut_cor = [abs(ide_cut(1:2:end-1)), sign(ide_cut(1:2:end-1)).*edge_bound_cut(id(1:2:end-1),1), sign(ide_cut(2:2:end)).*edge_bound_cut(id(2:2:end),1)];
    tri_cut_cor = [edge_bound_cut(id(1:2:end-1),2), edge_bound_cut(id(2:2:end),2)];
    ide_cut_cor(:,2:3) = (tri_cut_cor(:,1) == param.E2T(ide_cut_cor(:,1),1)).*ide_cut_cor(:,[2 3]) + (tri_cut_cor(:,2) == param.E2T(ide_cut_cor(:,1),1)).*ide_cut_cor(:,[3 2]);
    
    % Constraint rotations between each side of the cut
    R0 = eye(2,2);
    R1 = [0,-1; 1,0];
    R2 = R1^2;
    R3 = R1^3;
    R = [R0(:), R1(:), R2(:), R3(:)]';
    R = matrix_vector_multiplication(R(k21(ide_cut_cor(:,1)),:));
    I1 = sparse(1:size(ide_cut_cor,1), abs(ide_cut_cor(:,2)), sign(ide_cut_cor(:,2)), size(ide_cut_cor,1), SrcCut.ne);
    I2 = sparse(1:size(ide_cut_cor,1), abs(ide_cut_cor(:,3)), sign(ide_cut_cor(:,3)), size(ide_cut_cor,1), SrcCut.ne);
    Rot = blkdiag(I1*dec_cut.d0p,I1*dec_cut.d0p) - R*blkdiag(I2*dec_cut.d0p,I2*dec_cut.d0p);
    
    % Boundary edges are vertical or horizontal
    Align = sparse(0,2*SrcCut.nv);
    if (ifboundary && ~isempty(param.ide_bound)) || (ifhardedge && ~isempty(param.ide_hard))
        ide_fix_cut = find(ismember(abs(ide_cut_inv), param.ide_fix));
        tri_fix_cut = vec(SrcCut.E2T(ide_fix_cut,1:2));
        id = tri_fix_cut ~= 0;
        tri_fix_cut = tri_fix_cut(id);
        [~,ia] = ismember(param.tri_fix, tri_fix_cut);
        assert(length(ia) == length(tri_fix_cut))
        ide_fix_cut = [ide_fix_cut; ide_fix_cut];
        ide_fix_cut = ide_fix_cut(id);
        ide_fix_cut = ide_fix_cut(ia);
    
        dir_fix = round(abs(wrapToPi(2*ang(param.tri_fix))/pi) + 1);
        Align = sparse(1:length(ide_fix_cut), ide_fix_cut + (2 - dir_fix)*SrcCut.ne, 1, length(ide_fix_cut), 2*SrcCut.ne)*blkdiag(dec_cut.d0p, dec_cut.d0p);
    end
else
    Align = sparse(0,2*SrcCut.nv);
    Rot = sparse(0,2*SrcCut.nv);
end
