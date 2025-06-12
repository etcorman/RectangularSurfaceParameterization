function dec = dec_tri(Src)

nv = Src.nv;
ne = Src.ne;
nf = Src.nf;

l2 = [sum((Src.X(Src.T(:,2),:) - Src.X(Src.T(:,3),:)).^2, 2), ...
      sum((Src.X(Src.T(:,3),:) - Src.X(Src.T(:,1),:)).^2, 2), ...
      sum((Src.X(Src.T(:,1),:) - Src.X(Src.T(:,2),:)).^2, 2)];
half_area = l2.*abs(Src.cot_corner_angle)/4; % force positive area
half_area = (half_area(:,[2 3 1]) + half_area(:,[3 1 2]))/2;
half_area = repmat(Src.area/3, [1,3]);
vor_area = accumarray(Src.T(:), half_area(:));
cotweight = accumarray(abs(Src.T2E(:)), vec(Src.cot_corner_angle(:,[3 1 2]))/2);
assert(all(vor_area > 0), 'Negative vertex area.');

if any(cotweight < 1e-5)
    warning('Non Delaunay tet-mesh: risk of convergence issues!');
    cotweight = max(cotweight, 1e-5); % clamp to avoid problems
end

d0p = sparse([1:ne;1:ne]', Src.E2V, [ones(ne,1),-ones(ne,1)], ne, nv);
d1p = sparse([1:nf;1:nf;1:nf]', abs(Src.T2E), sign(Src.T2E), nf, ne);
assert(norm(d1p*d0p, 'fro') == 0, 'Assembling DEC: Orinetation problems');

star0p = sparse(1:nv, 1:nv, vor_area, nv, nv);
star1p = sparse(1:ne, 1:ne, cotweight, ne, ne);
star2p = sparse(1:nf, 1:nf, 1./Src.area, nf, nf);

d0d = d1p';
d1d = d0p';
assert(norm(d1d*d0d, 'fro') == 0, 'Assembling DEC: Orinetation problems');

star0d = sparse(1:nf, 1:nf, Src.area, nf, nf);
star1d = sparse(1:ne, 1:ne, 1./cotweight, ne, ne);
star2d = sparse(1:nv, 1:nv, 1./vor_area, nv, nv);

% output
dec.d0p = d0p;
dec.d1p = d1p;
dec.d0d = d0d;
dec.d1d = d1d;

dec.star0p = star0p;
dec.star1p = star1p;
dec.star2p = star2p;
dec.star0d = star0d;
dec.star1d = star1d;
dec.star2d = star2d;

% Triangle based operator
Reduction_tri = sparse(reshape((1:3*nf)', [nf,3]), Src.T, 1, 3*nf, nv);
deg_ed = accumarray(abs(Src.T2E(:)), 1);
I = abs(Src.T2E(:,[1 2 3])); 
J = reshape((1:3*nf),[nf,3]); 
S = sign(Src.T2E(:,[1 2 3]))./deg_ed(abs(Src.T2E));
d0p_tri = sparse([I, I], [J, J(:,[2 3 1])], [-S, S], ne, 3*nf);

star0p_tri = sparse(J, J, half_area, 3*nf, 3*nf);

W = d0p'*star1p*d0p;
W_tri = d0p_tri'*star1p*d0p_tri;

dec.W = (W + W')/2;
dec.d0p_tri = d0p_tri;
dec.star0p_tri = star0p_tri;
dec.W_tri = (W_tri + W_tri')/2;
dec.Reduction_tri = Reduction_tri;

end
