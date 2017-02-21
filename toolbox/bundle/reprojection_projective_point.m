function [x] = reprojection_projective_point( a, b )
% REPROJECTION_PROJECTIVE_POINT reprojects a point.
%
% x returns 2-by-1 vector

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

x_ = reshape(a,3,4) * [ b ; 1 ];
x = x_(1:2) / x_(3);

end
