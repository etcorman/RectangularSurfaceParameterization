function [omega,ang,sing] = trivial_connection(Src, param, dec, ifboundary, ifhardedge, sing, om_cycle, om_link)

% Default value for vertex indices
if ~exist('sing','var')
    if ~isempty(param.idx_bound)
        sing = zeros(Src.nv,1);
        sing(param.idx_bound) = round(4*param.Kt(param.idx_bound)/(2*pi))/4;
    else
        sing = zeros(Src.nv,1);
        id = randi(Src.nv, 4*(Src.nf - Src.ne + Src.nv), 1);
        sing(id) = 1/4;
    end
end
% Default value for non-contractible indices
if ~exist('om_cycle','var') || isempty(om_cycle)
    om_cycle = param.Icycle*param.para_trans;
    om_cycle = om_cycle - 2*pi*round(4*om_cycle/(2*pi))/4;
end
% Default value for connecting constraint path indices
if ~exist('om_link','var') || isempty(om_link)
    om_link = param.Ilink*param.para_trans;
    om_link = om_link - 2*pi*round(4*om_link/(2*pi))/4;
end

% Check Gauss-Bonnet constraints
if isempty(param.idx_bound)
    assert(norm(sum(sing) - (Src.nf - Src.ne + Src.nv)) < 1e-5, 'Singularities do not satisfy Gauss-Bonnet.');
end

% Build loop constraints
if ifboundary && ifhardedge
    s2 = [sing; round(2*param.K(Src.nv+1:end)/pi)/4];
    A = [param.d1d; sparse(1:length(param.ide_bound), param.ide_bound, 1, length(param.ide_bound), Src.ne); param.Ilink; param.Icycle];
    b = [param.K-2*pi*s2; zeros(length(param.ide_bound),1); om_link; om_cycle];
elseif ifboundary
    A = [dec.d1d; sparse(1:length(param.ide_bound), param.ide_bound, 1, length(param.ide_bound), Src.ne); param.Ilink; param.Icycle];
    b = [param.Kt-2*pi*sing; zeros(length(param.ide_bound),1); om_link; om_cycle];
elseif ifhardedge
    idx = setdiff((1:size(param.d1d,1))', param.idx_bound);
    s2 = [sing; round(2*param.K(Src.nv+1:end)/pi)/4];
    A = [param.d1d(idx,:); sparse(1:length(param.ide_hard), param.ide_hard, 1, length(param.ide_hard), Src.ne); param.Ilink; param.Icycle];
    b = [param.K(idx)-2*pi*s2(idx); zeros(length(param.ide_hard),1); om_link; om_cycle];
else
    A = [dec.d1d(param.idx_int,:); sparse(1:length(param.ide_bound), param.ide_bound, 1, length(param.ide_bound), Src.ne); param.Ilink; param.Icycle];
    b = [param.Kt(param.idx_int)-2*pi*sing(param.idx_int); zeros(length(param.ide_bound),1); om_link; om_cycle];
end

% Solve qudratic program
omega = quadprog(dec.star1d, zeros(Src.ne,1), [], [], A, b);

% Compute frame angle
ang = brush_frame_field(param, omega, param.tri_fix);

% Check that alignment constraints are satisfied
% assert(norm(wrapToPi(4*ang)) < 1e-5, 'Failed to prescribe alignment constraint.');

% Check that singularities indices are as prescribed
sing2 = (dec.d1d*(param.para_trans - omega) + param.Kt_invisible)/(2*pi);
assert(norm(sing(param.idx_int) - sing2(param.idx_int)) < 1e-5, 'Failed to prescribe singularities.');

sing_loop = wrapToPi(om_cycle - param.Icycle*param.para_trans)/(2*pi);
sing_loop2 = wrapToPi(param.Icycle*(omega - param.para_trans))/(2*pi);
assert(norm(sing_loop - sing_loop2) < 1e-5, 'Failing cycle constraints.');

sing_link = wrapToPi(om_cycle - param.Ilink*param.para_trans)/(2*pi);
sing_link2 = wrapToPi(param.Ilink*(omega - param.para_trans))/(2*pi);
assert(norm(sing_link - sing_link2) < 1e-5, 'Failed to prescribe constraints between feature curves.');

