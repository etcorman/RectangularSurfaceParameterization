function [omega,ang,sing,kappa,Curv] = compute_curvature_cross_field(Src, param, dec, smoothing_iter, alpha)
% Compute curvature aligned cross field
% - omega: field rotation
% -   ang: field angle
% -  sing: field cingularities
% - kappa: principal curvatures
% -  Curv: symmetric curvature tensor

%% Compute curvature tensor
% "SECOND FUNDAMENTAL MEASURE OF GEOMETRIC SETS AND LOCAL APPROXIMATION OF CURVATURES"
% David Cohen-Steiner & Jean-Marie Morvan
% Compute dihedral angle
comp_angle = @(u,v,n) atan2(dot(cross(u,v,2),n,2), dot(u,v,2));
edge = Src.X(Src.E2V(:,2),:) - Src.X(Src.E2V(:,1),:);
edge_length = sqrt(sum(edge.^2,2));
edge = edge./edge_length;
Eedge = [edge(:,1).*edge(:,1), edge(:,1).*edge(:,2), edge(:,1).*edge(:,3), ...
         edge(:,2).*edge(:,1), edge(:,2).*edge(:,2), edge(:,2).*edge(:,3), ...
         edge(:,3).*edge(:,1), edge(:,3).*edge(:,2), edge(:,3).*edge(:,3)];

ide_int = all(Src.E2T(:,1:2) ~= 0,2);
dihedral_angle = zeros(Src.ne,1);
dihedral_angle(ide_int) = Src.E2T(ide_int,4).*comp_angle(Src.normal(Src.E2T(ide_int,1),:), Src.normal(Src.E2T(ide_int,2),:), edge(ide_int,:));

% Compute curvature tensor
Curv = zeros(Src.nf,4);
J = [0,-1; 1, 0];
K = zeros(Src.nf,1);
for i = 1:Src.nf
    idt = i;
    idt = find(any(ismember(Src.T, Src.T(idt,:)),2));
    ide = unique(abs(Src.T2E(idt,:)));
    E = [param.e1r(i,:)', param.e2r(i,:)'];
    
    A = sum(dihedral_angle(ide).*edge_length(ide).*Eedge(ide,:))/Src.area(i);
    A = reshape(A, [3,3]);
    A = (A + A')/2;
    A = J'*E'*A*E*J;
    A = (A + A')/2;

    Curv(i,:) = A(:);
    K(i) = det(A);
end
Curv = Curv(:,[1 2 4]);

%% Extract principal directions
dir_max = zeros(Src.nf,1); % principal direction
kappa = zeros(Src.nf,2); % principal curvature
for i = 1:Src.nf
    A = [Curv(i,1), Curv(i,2); Curv(i,2), Curv(i,3)];

    [V,D] = eig(A);

    dir_max(i) = complex(V(1,1), V(2,1));
    kappa(i,:) = diag(D);
end

%% Frame field smoothing
z = (dir_max./abs(dir_max)).^4; % Init cross field from principal direction
z_fix = ones(length(param.tri_fix),1); % reference frame is algned with constraint edge by construction
z(param.tri_fix) = z_fix; % alignment constraints

if smoothing_iter > 0
    % Build connection Laplacian
    I = [param.ide_int,param.ide_int];
    J = param.E2T(param.ide_int,1:2);
    rot = param.para_trans(param.ide_int);
    S = [exp(1i*4*rot/2); -exp(-1i*4*rot/2)];
    d0d_cplx = sparse(I, J, S, Src.ne, Src.nf);
    Wcon = d0d_cplx'*dec.star1d*d0d_cplx;
    Wcon = (Wcon + Wcon')/2;

    % Screen smoothing
    w = (abs(K) + 1e-3); % Gaussian curvature weight
    M = spdiags(alpha*w.*Src.area, 0, Src.nf, Src.nf); % Modified mass matrix
    A = Wcon + M;
    for i = 1:smoothing_iter
        if ~isempty(param.tri_fix)
            z(param.tri_free) = A(param.tri_free,param.tri_free)\(M(param.tri_free,param.tri_free)*z(param.tri_free) - A(param.tri_free,param.tri_fix)*z_fix);
        else
            z = A\(M*z);
        end
        z = z./abs(z);
    end
end

%% Extract angles in reference basis
% Compute frames rotation
ang = angle(z)/4;
omega = wrapToPi(4*(dec.d0d*ang + param.para_trans))/4;
omega(param.ide_bound) = 0;
sing = (dec.d1d*(param.para_trans - omega) + param.Kt_invisible)/(2*pi);

% Brush frame field
ang = brush_frame_field(param, omega, param.tri_fix, ang(param.tri_fix));

%% Match curvature with closest frame direction
z_max = (dir_max./abs(dir_max)).^2;
z_min = 1i*z_max;
[~,id] = min(abs([z_max, z_min] - exp(2*1i*ang)), [], 2);
kappa = [kappa(:,1).*(id == 1) + kappa(:,2).*(id == 2), ...
         kappa(:,1).*(id == 2) + kappa(:,2).*(id == 1)];

% % Find nearest angle
% q = 6;
% ang_curv = angle(z_max)/2 + (-q:q)*pi/2;
% z_curv = exp(1i*ang_curv);
% [~,id] = min(abs(z_curv - exp(1i*ang)), [], 2);
% ang_curv = sum(ang_curv.*(id == (1:2*q+1)), 2);
% i = round((ang_curv - ang)/(2*pi));
% ang_curv = ang_curv - i*2*pi;
% assert(max(abs(ang_curv - ang))/(2*pi) < 0.25, 'Fail curvature alignment.');

