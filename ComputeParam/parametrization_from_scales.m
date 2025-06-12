function [Xp,mu] = parametrization_from_scales(Src, SrcCut, dec_cut, param, ang, omega, ut, vt, Align, Rot)

if (exist('Align','var') && ~isempty(Align)) || (exist('Rot','var') && ~isempty(Rot))
    ifseamless_const = true;
else
    ifseamless_const = false;
end
if ifseamless_const
    if ~exist('Align','var')
        Align = sparse(2*SrcCut.nv,0);
    end
    if ~exist('Rot','var')
        Rot = sparse(2*SrcCut.nv,0);
    end
end

% Integrated scales per edge
u1_int = ut(:,[1 2 3])/2 + ut(:,[2 3 1])/2;
u2_int = vt(:,[1 2 3])/2 + vt(:,[2 3 1])/2;

expu = [exp(u1_int + u2_int),      zeros(Src.nf,3), ...
             zeros(Src.nf,3), exp(u1_int - u2_int)];

% Average frame on edge
omega(param.ide_bound) = 0;
e1 = exp(1i*ang);
e1_edge = [e1, e1, e1];
e1_edge = (e1_edge + exp(-1i*omega(abs(Src.T2E)).*sign(Src.T2E)).*e1_edge)/2; 
e2_edge = 1i*e1_edge;

% Edge of cut mesh on local basis
edge = dec_cut.d0p*SrcCut.X;
edge1 = edge(abs(SrcCut.T2E(:,1)),:);
edge2 = edge(abs(SrcCut.T2E(:,2)),:);
edge3 = edge(abs(SrcCut.T2E(:,3)),:);
edge_tri_cut = [complex(dot(param.e1r, edge1, 2), dot(param.e2r, edge1, 2)), complex(dot(param.e1r, edge2, 2), dot(param.e2r, edge2, 2)), complex(dot(param.e1r, edge3, 2), dot(param.e2r, edge3, 2))];

% Deformed edges inside triangle (ie \mu_{ij}^k)
sigma1_tri = real(conj(e1_edge).*edge_tri_cut);
sigma2_tri = real(conj(e2_edge).*edge_tri_cut);
mu1_tri = expu(:,1:3).*sigma1_tri + expu(:,4:6)  .*sigma2_tri;
mu2_tri = expu(:,7:9).*sigma1_tri + expu(:,10:12).*sigma2_tri;

% New edge vector (ie \mu_{ij})
mu = zeros(SrcCut.ne,2);
mu(:,1) = accumarray(abs(SrcCut.T2E(:)), mu1_tri(:))./accumarray(abs(SrcCut.T2E(:)), 1);
mu(:,2) = accumarray(abs(SrcCut.T2E(:)), mu2_tri(:))./accumarray(abs(SrcCut.T2E(:)), 1);

% Integration with/out seamless constraints
W = dec_cut.W + 1e-5*dec_cut.star0p;
W = (W + W')/2;
div_dX = dec_cut.d1d*dec_cut.star1p*mu;

if ifseamless_const
    Xp = quadprog(blkdiag(W,W),-div_dX(:), [], [], [Align; Rot], zeros(size(Align,1)+size(Rot,1),1));
    Xp = reshape(Xp, [SrcCut.nv,2]);
else
    Xp = W\div_dX;
end
