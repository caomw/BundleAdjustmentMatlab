function [x] = reprojection_point( Kparam, a, b, num_variableK )
% REPROJECTION_POINT reprojects a point given camera and point parameters.
%
% x returns 2-by-1 vector

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

% build K
if num_variableK == 1
    Kparam(1) = a(7);
    Kparam(2) = a(7);
elseif num_variableK == 4
    Kparam = a(7:end);
end
K = calibration_matrix(Kparam);

% reproject
x_ = K * ( vl_rodr(a(1:3)) * b + a(4:6) );
x = x_(1:2) / x_(3);

end
