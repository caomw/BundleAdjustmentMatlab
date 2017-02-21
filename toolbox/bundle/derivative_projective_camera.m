function dfdx = derivative_projective_camera( x0, dx, f0, b )
% DERIVATIVE_PROJECTIVE_CAMERA computes the partial derivative w.r.t projective camera parameters.

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

% 1st order forward approximation
h = 1e-10;
dx = dx / norm(dx);
x1 = x0 + h * dx;
f1 = reprojection_projective_point(x1, b);
dfdx = (f1 - f0) / h;

end