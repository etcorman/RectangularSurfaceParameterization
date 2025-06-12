function [F,Jf,Hf] = oracle_integrability_condition(Src, param, dec, omega, ut, vt, ang, lambda, Reduction, ide_free)
% Compute the integrability condition per edge and its derivative

% OUTPUT:
% - F: integrability condition evaluated on edges defined by the list "ide_free"
% - Jf: Jacobian of F wrt u, v, theta
% - Hf: second order derivative of F needed for the Newton method

if ~exist('ide_free','var')
    ide_free = param.ide_free;
end

d0d = dec.d0d;
d0d(param.ide_hard,:) = 0;
d0d(param.ide_bound,:) = 0;

[O,Or,dO] = omega_from_scale(Src, param, dec, ut, vt, ang, Reduction);

F = O(ide_free,:)*[ut(:); vt(:)] - omega(ide_free);
Jf = [Or(ide_free,:), dO(ide_free,:) - d0d(ide_free,:)];

% Second order derivative
if nargout >= 3
    cot_ang = Src.cot_corner_angle;
    cos_2ff1 = cos(2*ang + 2*param.ang_basis(:,1));
    sin_2ff1 = sin(2*ang + 2*param.ang_basis(:,1));
    cos_2ff2 = cos(2*ang + 2*param.ang_basis(:,2));
    sin_2ff2 = sin(2*ang + 2*param.ang_basis(:,2));
    cos_2ff3 = cos(2*ang + 2*param.ang_basis(:,3));
    sin_2ff3 = sin(2*ang + 2*param.ang_basis(:,3));

    vt1  = vt(:,1);
    vt2  = vt(:,2);
    vt3  = vt(:,3);

    lambda_full = zeros(Src.ne,1);
    lambda_full(ide_free) = lambda(1:length(ide_free));
    le = sign(Src.T2E).*lambda_full(abs(Src.T2E));
    
    dthS = 2*le(:,1).*cot_ang(:,3).*(-sin_2ff1.*(cot_ang(:,2).*(vt3 - vt1) + cot_ang(:,1).*(vt3 - vt2)) - cos_2ff1.*(vt2 - vt1)) + ...
           2*le(:,2).*cot_ang(:,1).*(-sin_2ff2.*(cot_ang(:,3).*(vt1 - vt2) + cot_ang(:,2).*(vt1 - vt3)) - cos_2ff2.*(vt3 - vt2)) + ...
           2*le(:,3).*cot_ang(:,2).*(-sin_2ff3.*(cot_ang(:,1).*(vt2 - vt3) + cot_ang(:,3).*(vt2 - vt1)) - cos_2ff3.*(vt1 - vt3));
    D2_th = spdiags(dthS, 0, Src.nf, Src.nf);

    I = repmat((1:Src.nf)', [1,3]);
    J = reshape((1:3*Src.nf)', [Src.nf,3]);
    dvS1 = le(:,1).*cot_ang(:,3).*(-cos_2ff1.*cot_ang(:,2) + sin_2ff1) + ...
           le(:,2).*cot_ang(:,1).*cos_2ff2.*(cot_ang(:,3) + cot_ang(:,2)) + ...
           le(:,3).*cot_ang(:,2).*(-cos_2ff3.*cot_ang(:,3) - sin_2ff3);
    dvS2 = le(:,1).*cot_ang(:,3).*(-cos_2ff1.*cot_ang(:,1) - sin_2ff1) + ...
           le(:,2).*cot_ang(:,1).*(-cos_2ff2.*cot_ang(:,3) + sin_2ff2) + ...
           le(:,3).*cot_ang(:,2).*cos_2ff3.*(cot_ang(:,1) + cot_ang(:,3));
    dvS3 = le(:,1).*cot_ang(:,3).*cos_2ff1.*(cot_ang(:,2) + cot_ang(:,1)) + ...
           le(:,2).*cot_ang(:,1).*(-cos_2ff2.*cot_ang(:,2) - sin_2ff2) + ...
           le(:,3).*cot_ang(:,2).*(-cos_2ff3.*cot_ang(:,1) + sin_2ff3);
    D_vth = sparse(I, J, [dvS1, dvS2, dvS3], Src.nf, 3*Src.nf);
    assert(max(abs(D_vth*vt(:) - dO'*lambda_full)) < 1e-6, 'Second derivative of constraints is invalid.');

    Hf = [sparse(3*Src.nf,3*Src.nf), sparse(3*Src.nf,3*Src.nf), sparse(3*Src.nf,Src.nf);...
          sparse(3*Src.nf,3*Src.nf), sparse(3*Src.nf,3*Src.nf),                  D_vth';...
            sparse(Src.nf,3*Src.nf),                     D_vth,                  D2_th];
    Red = blkdiag(Reduction, speye(Src.nf,Src.nf));
    Hf = Red'*Hf*Red;
    Hf = (Hf + Hf')/2;
end