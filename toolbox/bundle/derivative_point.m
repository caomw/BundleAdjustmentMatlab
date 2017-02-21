function dfdx = derivative_point( x0, dx, f0, K, a, num_variableK )
% DERIVATIVE_POINT computes the partial derivative w.r.t point parameters.

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

% 1st order forward approximation
h = 1e-10;
dx = dx / norm(dx);
x1 = x0 + h * dx;
f1 = reprojection_point(K, a, x1, num_variableK);
dfdx = (f1 - f0) / h;

end