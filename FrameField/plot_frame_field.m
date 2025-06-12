function fig = plot_frame_field(fig, Src, param, ang, col)

if isempty(fig)
    fig = figure;
end

figure(fig);
trisurf(Src.T, Src.X(:,1), Src.X(:,2), Src.X(:,3), col, 'edgecolor', 'none');
if size(ang,1) == size(Src.T,1)
    e1 = exp(1i*ang);
    e2 = 1i*e1;
    E1 = real(e1).*param.e1r + imag(e1).*param.e2r;
    E2 = real(e2).*param.e1r + imag(e2).*param.e2r;
    
    bar = (Src.X(Src.T(:,1),:) + Src.X(Src.T(:,2),:) + Src.X(Src.T(:,3),:))/3;
    
    hold on;
    quiver3(bar(:,1), bar(:,2), bar(:,3), E1(:,1), E1(:,2), E1(:,3), 0.5, 'r', 'LineWidth',1);
    quiver3(bar(:,1), bar(:,2), bar(:,3), E2(:,1), E2(:,2), E2(:,3), 0.5, 'g', 'LineWidth',1);
    hold off;
end
axis equal; view(0,90); % caxis([-0.5 0.5])

if size(col,1) == Src.nv
    shading interp;
end