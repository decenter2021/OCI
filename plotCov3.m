function [x,y,z] = plotCov3(Y)
    [V, D] = eig(Y);
    [u, v] = meshgrid(linspace(0, 2*pi, 100), linspace(0, pi, 50));
    x_sphere = cos(u) .* sin(v);
    y_sphere = sin(u) .* sin(v);
    z_sphere = cos(v);
    ellipsoid_points = V * (D^(-1/2)) * [x_sphere(:)'; y_sphere(:)'; z_sphere(:)'];
    x = reshape(ellipsoid_points(1, :), size(x_sphere));
    y = reshape(ellipsoid_points(2, :), size(y_sphere));
    z = reshape(ellipsoid_points(3, :), size(z_sphere));
end