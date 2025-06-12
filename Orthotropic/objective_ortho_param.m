function [fct,H,df] = objective_ortho_param(energy_type, weight, Src, dec, param, angn, ut, vt, Reduction)

Wt = dec.W_tri;
AeT = dec.star0p_tri;
Iso = blkdiag(AeT, AeT);
Conf = blkdiag(0*AeT, AeT);
Ar = blkdiag(AeT, 0*AeT);
integral_tri = ones(3,1)/3;

H = sparse(2*Src.nv+Src.nf,2*Src.nv+Src.nf);
if strcmp(energy_type, 'distortion')
    % H_u = (1 - weight.w_conf_ar)*Ar + weight.w_conf_ar*Conf;
    H_u = (1 - weight.w_conf_ar)*AeT;
    H_v =       weight.w_conf_ar*AeT;
    H_ang = sparse(Src.nf, Src.nf);

    H_uv = blkdiag(H_u, H_v);
    df = [Reduction'*[H_u*ut(:); H_v*vt(:)]; H_ang*angn];
    fct = [ut(:); vt(:)]'*[H_u*ut(:); H_v*vt(:)];
elseif strcmp(energy_type, 'chebyshev')
    % Chebyshev Net: |J\[1;1]| = sqrt(2) (unit diagonal on the shape)
    % linearization of log(exp(-2*a-2*b)/2+exp(-2*a+2*b)/2) = 0
    % err_diag = log(exp(2*ut*[-integral_tri;-integral_tri])/2 + exp(2*ut*[-integral_tri; integral_tri])/2);
    err_diag = log(exp(- 2*ut*integral_tri - 2*vt*integral_tri)/2 + exp(- 2*ut*integral_tri + 2*vt*integral_tri)/2);

    % Jacobian
    da = -2*ones(Src.nf,1)/3;
    db = (2/3)*(exp(4*vt*integral_tri) - 1)./(exp(4*vt*integral_tri) + 1);

    I = repmat((1:Src.nf)', [1,3]);
    J = reshape((1:3*Src.nf)', [Src.nf,3]);
    Da = sparse(I, J, [da,da,da], Src.nf, 3*Src.nf);
    Db = sparse(I, J, [db,db,db], Src.nf, 3*Src.nf);
    
    % Second order
    daa = zeros(Src.nf,1);
    dab = zeros(Src.nf,1);
    dbb = (16/9)*exp(4*vt*integral_tri)./(exp(4*vt*integral_tri) + 1).^2;

    Haa = sparse([I,Src.nf+I,2*Src.nf+I], [J,J,J], repmat(daa.*err_diag.*Src.area, [1,9]), 3*Src.nf, 3*Src.nf);
    Hab = sparse([I,Src.nf+I,2*Src.nf+I], [J,J,J], repmat(dab.*err_diag.*Src.area, [1,9]), 3*Src.nf, 3*Src.nf);
    Hbb = sparse([I,Src.nf+I,2*Src.nf+I], [J,J,J], repmat(dbb.*err_diag.*Src.area, [1,9]), 3*Src.nf, 3*Src.nf);

    H_uv = [Da, Db]'*dec.star0d*[Da, Db] + [Haa, Hab; Hab', Hbb];
    H_ang = sparse(Src.nf, Src.nf);

    df = [Reduction'*([Da, Db]'*dec.star0d*err_diag); zeros(Src.nf,1)];
    fct = err_diag'*dec.star0d*err_diag;
elseif strcmp(energy_type, 'alignment')
    % Stay close to given direction
    ang_dir = weight.ang_dir;

    % Aspect ratio constraint
    log_aspect_ratio = log(weight.aspect_ratio)/2;
    diff = [ut(:); vt(:)] - repmat(log_aspect_ratio, [6,1]);
    H_uv = weight.w_ratio*Conf;
    df_uv = Reduction'*H_uv*diff;
    fct = diff'*H_uv*diff;

    % Angle constraint
    H_ang = weight.w_ang*dec.star0d;
    df_ang = H_ang*(angn - ang_dir);
    fct = fct + (angn - ang_dir)'*(H_ang*(angn - ang_dir));

    df = [df_uv; df_ang];
else
    error('This energy does not exist.');
end

if isfield(weight,'w_gradv') && (weight.w_gradv > 0)
    H_uv = H_uv + weight.w_gradv*blkdiag(0*Wt, Wt);
    fct = fct + weight.w_gradv*vt(:)'*Wt*vt(:);
    df = df + [weight.w_gradv*Reduction'*blkdiag(0*Wt, Wt)*[ut(:); vt(:)]; zeros(Src.nf,1)];
end

H = blkdiag(Reduction'*H_uv*Reduction, H_ang);
H = (H + H')/2;

