function [O,Or,dO] = omega_from_scale(Src, param, dec, ut, vt, ang, Reduction)
% Build frame rotation from scale factors ut and vt defined at corners

% OUTPUT:
% - O: (E x 6F sparse matrix) O*[ut(:); vt(:)] is the frame rotation omega
% - Or: (E x 2V sparse matrix) linear map from scale factors to omega in
% reduced variables
% - dO: (E x F sparse matrix) derivative of omega wrt ang

cot_ang = Src.cot_corner_angle;

cos_2ff1 = cos(2*ang + 2*param.ang_basis(:,1));
sin_2ff1 = sin(2*ang + 2*param.ang_basis(:,1));
cos_2ff2 = cos(2*ang + 2*param.ang_basis(:,2));
sin_2ff2 = sin(2*ang + 2*param.ang_basis(:,2));
cos_2ff3 = cos(2*ang + 2*param.ang_basis(:,3));
sin_2ff3 = sin(2*ang + 2*param.ang_basis(:,3));

I = Src.T2E(:,[1 1 1 2 2 2 3 3 3]);
J = repmat(reshape((1:3*Src.nf)', [Src.nf,3]), [1,3]);
S = 0.5*sign(I).*[cot_ang(:,3).*[-cos_2ff1 - sin_2ff1.*cot_ang(:,2)      , cos_2ff1 - sin_2ff1.*cot_ang(:,1)      , sin_2ff1.*(cot_ang(:,1) + cot_ang(:,2))], ...
                  cot_ang(:,1).*[ sin_2ff2.*(cot_ang(:,2) + cot_ang(:,3)),-cos_2ff2 - sin_2ff2.*cot_ang(:,3)      , cos_2ff2 - sin_2ff2.*cot_ang(:,2)      ], ...
                  cot_ang(:,2).*[ cos_2ff3 - sin_2ff3.*cot_ang(:,3)      , sin_2ff3.*(cot_ang(:,3) + cot_ang(:,1)),-cos_2ff3 - sin_2ff3.*cot_ang(:,1)      ]];
Dv_tri = sparse(abs(I), J, S, Src.ne, 3*Src.nf);

O = [-dec.star1p*dec.d0p_tri, Dv_tri];

if exist('Reduction','var') && ~isempty(Reduction)
    Or = O*Reduction;
end

% Derivative wrt frame angle
if nargout >= 3
    vt1  = vt(:,1);
    vt2  = vt(:,2);
    vt3  = vt(:,3);
    
    I = Src.T2E;
    J = repmat((1:Src.nf)', [1,3]);
    S = sign(I).*[cot_ang(:,3).*(cos_2ff1.*(cot_ang(:,2).*(vt3 - vt1) + cot_ang(:,1).*(vt3 - vt2)) - sin_2ff1.*(vt2 - vt1)), ...
                  cot_ang(:,1).*(cos_2ff2.*(cot_ang(:,3).*(vt1 - vt2) + cot_ang(:,2).*(vt1 - vt3)) - sin_2ff2.*(vt3 - vt2)), ...
                  cot_ang(:,2).*(cos_2ff3.*(cot_ang(:,1).*(vt2 - vt3) + cot_ang(:,3).*(vt2 - vt1)) - sin_2ff3.*(vt1 - vt3))];
    dO = sparse(abs(I), J, S, Src.ne, Src.nf);
end
