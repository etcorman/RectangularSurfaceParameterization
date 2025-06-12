% Rectangle Surface Parametrization
% Corman and Crane, 2025

clear all; close all; clc;
addpath(genpath(pwd));

path_save = 'Results/';
path_data = 'Mesh/';

%% Options
% Select mesh
mesh_name = 'B36'; % SquareMyles B36 pig

% Frame field options
frame_field_type = 'smooth'; % 'curvature' 'smooth' 'trivial'
ifhardedge   = true;     % Constrain hard edges
ifboundary   = true;     % Aligned with boundary edges
ifseamless_const = true; % Force seamlessness when building parametrization
ifquantization = true;   % Compute integer seamless parametrization

% Select objective function
energy_type = 'distortion';  %  'distortion' 'chebyshev' 'alignment'

% Energy weights
if strcmp(energy_type, 'distortion')
    weight.w_conf_ar = 0.5; % 0: area-preserving -- 0.5 : isometry -- 1: conformal
end
if strcmp(energy_type, 'alignment')
    weight.w_ang   = 1; % weight direction energy
    weight.w_ratio = 1; % weight aspect ratio energy
end
weight.w_gradv = 1e-2;  % weight regularization on v

%% Load mesh
[X,T] = readOBJ([path_data, mesh_name, '.obj']);

% Rescale: area equals one
area_tot = sum(sqrt(sum(cross(X(T(:,1),:) - X(T(:,2),:), X(T(:,1),:) - X(T(:,3),:),2).^2,2)))/2;
X = X/sqrt(area_tot);

% Preprocess geometry
Src = MeshInfo(X, T);   % Contains connectivity and basic info
dec = dec_tri(Src);     % Compute DEC matrices
[param,Src,dec] = preprocess_ortho_param(Src, dec, ifboundary, ifhardedge, 40); % Constraints and remesh

col = zeros(Src.nv,1); col(Src.E2V(param.ide_fix,:)) = 1;
figure;
trisurf(Src.T, Src.X(:,1), Src.X(:,2), Src.X(:,3), col, 'facecolor','interp', 'edgecolor','k');
axis equal; title('Constraint')

%% Compute initial cross field
if strcmp(frame_field_type, 'curvature') % Follow curvature
    % Cross field from smoothed curvature
    [omega,ang,sing,kappa,Curv] = compute_curvature_cross_field(Src, param, dec, 30, 1e-1);

    % Target aspect ratio
    weight.aspect_ratio = ((abs(kappa(:,1)) + 1e-5)./(abs(kappa(:,2)) + 1e-5));
    t = exp(5); 
    weight.aspect_ratio = max(min(weight.aspect_ratio, t), 1/t);

    % Target direction
    weight.ang_dir = ang;

    figure; 
    trisurf(Src.T, Src.X(:,1), Src.X(:,2), Src.X(:,3), abs(log(weight.aspect_ratio))); 
    axis equal; title('Aspect ratio');
elseif strcmp(frame_field_type, 'smooth') % Compute smooth cross field
    [omega,ang,sing] = compute_face_cross_field(Src, param, dec, 10);
elseif strcmp(frame_field_type, 'trivial') % Trivial connection (here with no inner singularities)
    % Vertex singularity index
    sing = zeros(Src.nv,1);
    sing(param.idx_bound) = round(2*param.K(param.idx_bound)/pi)/4;

    % Singularity index of non-contractible cycles
    om_cycle = param.Icycle*param.para_trans;
    om_cycle = om_cycle - 2*pi*round(4*om_cycle/(2*pi))/4;

    % Singularity index between disconnected constraints
    om_link = param.Ilink*param.para_trans;
    om_link = om_link - 2*pi*round(4*om_link/(2*pi))/4;

    % Trivial connection
    [omega,ang,sing] = trivial_connection(Src, param, dec, ifboundary, ifhardedge, sing);
else
    error('Cross field option unavailable.')
end

% Plot frame field
plot_frame_field(1, Src, param, ang, sing);
title('Init frame field'); 

%% Compute cross field jumps and build reduction matrix for v
% Order triangles around each vertex to compute sign bits
[Edge_jump,v2t,base_tri] = reduce_corner_var_2d(Src); 

% Compute reduction matrix from corner variables to vertex variable
[k21,Reduction] = reduction_from_ff2d(Src, param, ang, omega, Edge_jump, v2t); 

%% Optimize integrability condition
% Solver option
itmax = 200;
ifplot = false;

% Init u
u = zeros(Src.nv,1);
v = zeros(Src.nv,1);

% Optimize with Newton method
[u,v,ut,vt,om,angn,flag] = optimize_RSP(omega, ang, u, v, Src, param, dec, Reduction, energy_type, weight, ifplot, itmax);

%% Compute parametrization
% Cut geometry
[SrcCut,dec_cut,Align,Rot] = mesh_to_disk_seamless(Src, param, angn, sing, k21, ifseamless_const, ifboundary, ifhardedge);

% Reconstruct parametrization
[Xp,dX] = parametrization_from_scales(Src, SrcCut, dec_cut, param, angn, om, ut, vt, Align, Rot);

%% Plot stuff
disto = extract_scale_from_param(Xp, Src.X, Src.T, param, SrcCut.T, angn);
curl_dX = sqrt(sum((dec_cut.d1p*dX).^2,2))./Src.area;

figure;
subplot(2,2,1);
trisurf(SrcCut.T, SrcCut.X(:,1), SrcCut.X(:,2), SrcCut.X(:,3), log10(curl_dX), 'edgecolor', 'none');
axis equal; view(0,-90); colorbar;
title('Integrability');
subplot(2,2,2);
trisurf(SrcCut.T, Xp(:,1), Xp(:,2), 0*Xp(:,2));
axis equal; view(0,-90);
title('Param');
subplot(2,2,3);
trisurf(SrcCut.T, SrcCut.X(:,1), SrcCut.X(:,2), SrcCut.X(:,3), log10(disto.area), 'edgecolor', 'none');
axis equal; view(0,-90); colorbar;
title('log area');
subplot(2,2,4);
trisurf(SrcCut.T, SrcCut.X(:,1), SrcCut.X(:,2), SrcCut.X(:,3), abs(log10(disto.conf)), 'edgecolor', 'none');
axis equal; view(0,-90); colorbar;
title('log conformal');

% Plot singularities
col = zeros(Src.nf,1); col(disto.detJ <= 0) = 1;
id_sing_p = sing > 1/8;
id_sing_m = sing <-1/8;
figure;
hold all;
trisurf(Src.T, Src.X(:,1), Src.X(:,2), Src.X(:,3), col);
scatter3(Src.X(id_sing_p,1), Src.X(id_sing_p,2), Src.X(id_sing_p,3), 100, 'r', 'filled');
scatter3(Src.X(id_sing_m,1), Src.X(id_sing_m,2), Src.X(id_sing_m,3), 100, 'b', 'filled');
hold off;
axis equal; view(0,-90); 
title([num2str(sum(id_sing_p)+sum(id_sing_m)), ' singus']);

%% Save mesh
if contains(energy_type, 'cheby')
    r = [1,1;-1,1]*(sqrt(2)/2);
    UV = Xp*r;
else
    UV = Xp;
end

if ifquantization && ~isempty(param.ide_bound) && ~ifboundary
    warning('Boundary alignment is disactivated.');
    warning('The Quantization step does not support free boundaries.');
end

if ifquantization && ~ifseamless_const
    warning('The seamlessness constraint is disactivated.');
    warning('The Quantization needs an exact seamless map as input.');
end

if ifquantization && any(disto.detJ < 0)
    warning('The parametrization is not locally injective.');
    warning('The Quantization step will fail.');
end
save_param(ifquantization, path_save, mesh_name, sqrt(area_tot)*Src.X, Src.T, UV, SrcCut.T, sing, Src.E2V(param.ide_hard,:));
