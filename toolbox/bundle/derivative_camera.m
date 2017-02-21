function dfdx = derivative_camera( x0, dx, f0, K, b, num_variableK )
% DERIVATIVE_CAMERA computes the partial derivative w.r.t camera parameters.

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

% 1st order forward approximation
h = 1e-10;
dx = dx / norm(dx);
x1 = x0 + h * dx;
f1 = reprojection_point(K, x1, b, num_variableK);
dfdx = (f1 - f0) / h;

end