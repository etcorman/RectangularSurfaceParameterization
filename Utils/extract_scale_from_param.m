function [disto,ut,theta,u_tri] = extract_scale_from_param(Xp, X, T, param, T_cut, ang)

nf = size(T,1);

% Compute distortion
u_tri = zeros(nf,2);
disto.ang_param = zeros(nf,1);
disto.detJ = zeros(nf,1);
disto.area = zeros(nf,1);
disto.conf = zeros(nf,1);
disto.orth = zeros(nf,1);
disto.cheb = zeros(nf,1);
disto.cheb_ang = zeros(nf,1);
for i = 1:nf
    Jac = inv([param.e1r(i,:); param.e2r(i,:)]*[X(T(i,1),:) - X(T(i,2),:); X(T(i,1),:) - X(T(i,3),:)]');
    Jac = [Xp(T_cut(i,1),1:2) - Xp(T_cut(i,2),1:2); Xp(T_cut(i,1),1:2) - Xp(T_cut(i,3),1:2)]'*Jac;

    [U,S,V] = svd(Jac);
    s = diag(S);

    Q = U*V';
    if det(Q) < 0
        d = [1,1]; [~,id] = min(s); d(id) =-1;
        Q = U*diag(d)*V';
    end
    disto.ang_param(i) = atan2(Q(1,2), Q(1,1));

    % u_tri(i,1) = log(s(1)*s(2))/2;
    % u_tri(i,2) = log(s(1)/s(2))/2;
    d = diag(U*S*U');
    u_tri(i,1) = log(d(1)*d(2))/2;
    u_tri(i,2) = log(d(1)/d(2))/2;
    if exist('ang','var') && abs(Q(1,:)*[cos(ang(i)); sin(ang(i))]) < 0.5
        u_tri(i,2) =-u_tri(i,2);
    end

    disto.area(i) = s(1)*s(2);
    disto.conf(i) = s(1)/s(2);
    disto.detJ(i) = det(Jac);
    disto.orth(i) = acos((Jac(1,:)*Jac(2,:)')/(norm(Jac(1,:))*norm(Jac(2,:))));

    u = (Jac\[1;1])/sqrt(2); v = (Jac\[-1;1])/sqrt(2);
    disto.cheb(i) = ((norm(u) - 1)^2 + (norm(v) - 1)^2);
    disto.cheb_ang(i) = acos((u'*v)/(norm(u)*norm(v)))*180/pi;
end
if any(disto.detJ <= 0)
    warning([num2str(sum(disto.detJ <= 0)), ' negative determinant.']);
end

if nargout > 1
    theta = angle(exp(4*1i*(disto.ang_param - ang)))/4;
    % theta(param.tri_fix) = 0;
    
    % Average on vertices
    u = accumarray(T(:), [u_tri(:,1); u_tri(:,1); u_tri(:,1)])./accumarray(T(:), 1);
    v = accumarray(T_cut(:), [u_tri(:,2); u_tri(:,2); u_tri(:,2)])./accumarray(T_cut(:), 1);
    ut = [u(T), v(T)];
end