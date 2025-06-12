function [omega,ang,sing] = compute_face_cross_field(Src, param, dec, smoothing_iter)
% Compute amooth cross field
% - omega: field rotation
% -   ang: field angle
% -  sing: field cingularities

power = 4;

%% Change connection
% If hardedge or boundary take care of acute angles
% cf "Frame Fields for CAD models", Advances in Visual Computing, 2021
% https://inria.hal.science/hal-03537852/
if ~isempty(param.tri_fix)
    ifcadff = true;
else
    ifcadff = false;
end

% Compute new connection transforming acute angles into right angle
if ifcadff
    tol = pi/16;

    % Find indices of acute corner angle
    id = param.K(param.idx_fix_plus) > pi/2;
    idx = param.Vp2V(param.idx_fix_plus(id));
    idxp = param.Vp2V(param.idx_fix_plus);

    % Deform Gaussian curvature
    d = (Src.X(idx,1)' - Src.X(idxp,1)).^2 + (Src.X(idx,2)' - Src.X(idxp,2)).^2 + (Src.X(idx,3)' - Src.X(idxp,3)).^2;
    K_new = param.K(param.idx_fix_plus);
    K_new((K_new >-tol) & (K_new < tol) & any(d < 1e-3*repmat(max(d,[],1), [size(d,1),1]),2)) = 0;
    K_new(id) = pi/2;

    % Find new connection
    H = dec.star1d + 1e-3*(dec.d1d'*dec.star2d*dec.d1d); H = (H' + H)/2;
    A = sparse(1:length(param.ide_fix), param.ide_fix, 0, length(param.ide_fix), Src.ne);
    b = zeros(length(param.ide_fix),1);
    omega_cadff = quadprog(H,[], [], [], [A; param.d1d(param.idx_fix_plus,:)], [b; param.K(param.idx_fix_plus) - K_new], [], [], []);
end

%% Compute cross field
% Build connection Laplacian
I = [param.ide_int,param.ide_int];
J = param.E2T(param.ide_int,1:2);
if ifcadff
    rot = param.para_trans(param.ide_int) - omega_cadff(param.ide_int);
else
    rot = param.para_trans(param.ide_int);
end
S = [exp(1i*power*rot/2); -exp(-1i*power*rot/2)];
d0d_cplx = sparse(I, J, S, Src.ne, Src.nf);
Wcon = d0d_cplx'*dec.star1d*d0d_cplx;
Wcon = (Wcon + Wcon')/2;

% Set constraints
tri_fix = param.tri_fix;
z_fix = ones(length(tri_fix),1); % reference frame is algned with constraint edge by construction
tri_free = setdiff((1:Src.nf)', tri_fix);

% Compute initial cross field
z = zeros(Src.nf,1);
z(tri_fix) = z_fix;
if ~isempty(tri_fix) % If boundaries: solve Poisson problem
    z(tri_free) =-Wcon(tri_free,tri_free)\(Wcon(tri_free,tri_fix)*z_fix);
    D = eigs(Wcon, dec.star0d, 5, 'sm');
    dt = 20*real(D(2));
else % If no boundaries: compute smallest eigenvector
    [P,D] = eigs(Wcon, dec.star0d, 5, 'sm');
    z = P(:,1);
    dt = 20*real(D(1,1));
end
z = z./abs(z);

% Smoothing by heat flow
A = Wcon + dt*dec.star0d;
for i = 1:smoothing_iter
    if ~isempty(tri_fix)
        z(tri_free) = A(tri_free,tri_free)\(dt*dec.star0d(tri_free,tri_free)*z(tri_free) - A(tri_free,tri_fix)*z_fix);
    else
        z = A\(dt*dec.star0d*z);
    end
    z = z./abs(z);
end
assert(all(~isnan(z)), 'NaN vector field.');

%% Extract angles in reference basis
% Compute rotation
ang = angle(z)/power;
if ifcadff
    omega = wrapToPi(power*(dec.d0d*ang + param.para_trans - omega_cadff))/power + omega_cadff;
else
    omega = wrapToPi(power*(dec.d0d*ang + param.para_trans))/power;
end
omega(param.ide_bound) = 0;
omega(param.ide_fix) = 0;

% Compute singularities
sing = (dec.d1d*(param.para_trans - omega) + param.Kt_invisible)/(2*pi);

% Brush cross field
ang = brush_frame_field(param, omega, tri_fix, ang(tri_fix));
