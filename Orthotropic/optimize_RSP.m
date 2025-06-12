function [u,v,ut,vt,om,angn,flag,it] = optimize_RSP(omega, ang, u, v, Src, param, dec, Reduction, energy_type, weight, if_plot, itmax, A_const, b_const)
% Newton solver minimizing the enrgy defined in "energy_type" subject to
% the satisfaction of the integrability constraint

% flag : 1 reached max iter, 0 convergence, -1 linesearch failed

% Default max iter
if ~exist('itmax', 'var')
    itmax = 300;
end

% Default linear constraints
if ~exist('A_const', 'var') || (size(A_const,2) ~= size(Reduction,2) + Src.nf)
    dof = size(Reduction,2);
    A_const = sparse(1:length(param.tri_fix), dof+param.tri_fix, 1, length(param.tri_fix), dof+Src.nf);
    b_const = zeros(length(param.tri_fix),1); % Assume that the cross field already respect bnd conditions
end

% Line-search parameters
rho = 0.9;
beta = 1;

% List of edges respecting integrability condition
ide_free = setdiff((1:Src.ne)', param.ide_fix);

% Initilization of optimization variables
om = omega;
angn = ang;
ut = reshape(Reduction*[u; v], [Src.nf,6]);
vt = ut(:,4:6);
ut = ut(:,1:3);

n_dual = length(ide_free) + size(A_const,1);
lambda = zeros(n_dual,1);

%
d0d = dec.d0d;
d0d(param.ide_hard,:) = 0;
d0d(param.ide_bound,:) = 0;

% Keep track of objective evolution
flag = 1;
fct = zeros(itmax+1,1);
fct_const = zeros(itmax+1,2);
fct_grad_norm = zeros(itmax+1,2);
err = zeros(itmax+1,1);
for it = 1:itmax
    disp(['-- Iteration ', num2str(it)]);
    weight.om = om;

    % Compute derivate of constraints
    [F,Jf,Hf] = oracle_integrability_condition(Src, param, dec, om, ut, vt, angn, lambda, Reduction, ide_free);
    F = [F; A_const*[u; v; zeros(Src.nf,1)] - b_const];
    Jf = [Jf; A_const];

    fct_const(it,1) = norm(F);
    fct_const(it,2) = max(abs(F));

    % Compute derivate of objective functions
    [fct(it),Hfct,dfct] = objective_ortho_param(energy_type, weight, Src, dec, param, angn, ut, vt, Reduction);
    fct_grad_norm(it,1) = norm(dfct + Jf'*lambda);
    fct_grad_norm(it,2) = max(abs(dfct + Jf'*lambda));
    err(it) = sqrt(fct_const(it,1)^2 + fct_grad_norm(it,1)^2);

    % Newton on KKT conditions
    A = [Hfct + Hf, Jf'; Jf, sparse(n_dual,n_dual)]; A = (A + A')/2;
    b =-[dfct + Jf'*lambda; F];

    % Solve quadratic problem
    H = blkdiag(dec.star0p, dec.star0p, dec.star0d, dec.star1p(ide_free,ide_free), dec.star0d(param.tri_fix,param.tri_fix));
    H = (H + H')/2;
    f = zeros(size(H,1),1);

    x = A\b;
    try
        assert(norm(A*x - b) < 1e-5, 'Optimization failed.');
    catch
        x = quadprog(H, f, [], [], A, b);
        assert(norm(A*x - b) < 1e-5, 'Optimization failed.');
    end

    % Line-search
    run = true;
    t = 1;
    while run
        % Step size
        t = min(1, beta/err(it));

        % Update variables
        u_new = u + t*x(1:Src.nv);
        v_new = v + t*x(Src.nv+1:2*Src.nv);
        ut_new = reshape(Reduction*[u_new; v_new], [Src.nf,6]);
        vt_new = ut_new(:,4:6);
        ut_new = ut_new(:,1:3);
        alp_new = x(2*Src.nv+1:2*Src.nv+Src.nf);
        angn_new = angn + t*alp_new;
        om_new = om + t*d0d*alp_new;
        lambda_new = lambda + t*x(2*Src.nv+Src.nf+1:end);

        % Update constraints
        [F,Jf] = oracle_integrability_condition(Src, param, dec, om_new, ut_new, vt_new, angn_new, lambda_new, Reduction, ide_free);
        F = [F; A_const*[u_new; v_new; zeros(Src.nf,1)] - b_const];
        Jf = [Jf; A_const];

        % Update objective
        [fct(it+1),~,dfct] = objective_ortho_param(energy_type, weight, Src, dec, param, angn_new, ut_new, vt_new, Reduction);

        fct_const(it+1,1) = norm(F);
        fct_const(it+1,2) = max(abs(F));
        fct_grad_norm(it+1,1) = norm(dfct + Jf'*lambda_new);
        fct_grad_norm(it+1,2) = max(abs(dfct + Jf'*lambda_new));
        err(it+1) = sqrt(fct_const(it+1,1)^2 + fct_grad_norm(it+1,1)^2);

        % Check if the search is over
        if t == 1
            if_end_search = err(it+1) < err(it)^2/(2*beta);
        else
            if_end_search = err(it+1) < err(it) - beta/2;
        end

        if if_end_search
            % Found new variables
            run = false;
            u = u_new;
            ut = ut_new;
            v = v_new;
            vt = vt_new;
            angn = angn_new;
            alp = t*alp_new;
            om = om_new;
            lambda = lambda_new;

            beta = beta/rho;
        else
            % Keep the search going by reducing the step size
            beta = beta*rho;
        end

        if t < 1e-12
            warning('Linesearch failed.');
            break;
        end
    end
    if t < 1e-12
        flag =-1;
        break;
    end

    % Display optimization energies
    disp(['Total error : ', num2str(err(it+1,1)), ' -- Objective : ', num2str(fct(it+1,1))]);
    disp(['Grad norm : ', num2str(fct_grad_norm(it+1,1)), ' -- Max : ', num2str(fct_grad_norm(it+1,2))]);
    disp(['Integrability : ', num2str(fct_const(it+1,1)), ' -- Max : ', num2str(fct_const(it+1,2))]);
    err_ang = abs(alp)*180/pi;
    disp(['Max frame field angle change : ', num2str(max(abs(err_ang)))]);

    % Show new frame field
    if if_plot
        plot_frame_field(1, Src, param, angn, err_ang);
        title(['New frame field ', num2str(it)]); colorbar;
    end
    
    % Check that boundary constraints still hold
    if ~isempty(param.tri_fix)
        err_ang_bound = (180/pi)*wrapToPi(4*angn(param.tri_fix) - 4*ang(param.tri_fix))/4;
        assert(max(abs(err_ang_bound)) < 1e-3, 'Boundary constraints not respected.');
    end

    % Stop optimization when converged
    if (err(it+1) < 1e-5) && (max(abs(err_ang)) < 1e-3)
        flag = 0;
        break;
    end
    % pause
end
disp('-- end loop --');

% Show convergence plots
if if_plot
    figure;
    semilogy([fct_grad_norm(1:it+1,1), fct_grad_norm(1:it+1,2)], 'linewidth', 2);
    legend({'Grad norm L^2' 'Grad norm L^\infty'}, 'fontsize', 14);
    title('Grad norm');
    grid on;
    
    figure;
    semilogy([fct_const(1:it+1,1), fct_const(1:it+1,2)], 'linewidth', 2);
    legend({'Integrability L^2' 'Integrability L^\infty'}, 'fontsize', 14);
    title('Integrability');
    grid on;
    
    figure;
    plot(fct(1:it+1), 'linewidth', 2);
    title('Objective');
end