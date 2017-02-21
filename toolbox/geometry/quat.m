function R = quat( q )
% QUAT converts a quaternion to a rotation matrix
%
% Input:
%       q - 4x1 quaternion vector
%
% Output:
%       R - 3x3 rotation matrix

% reference: Conversion Quaternion to Matrix - Martin Baker

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

q2 = q .* q;
invs = 1 / sum(q2);

R = zeros(3,3);

R(1,1) = (q2(1) + q2(2) - q2(3) - q2(4)) * invs;
R(2,2) = (q2(1) - q2(2) + q2(3) - q2(4)) * invs;
R(3,3) = (q2(1) - q2(2) - q2(3) + q2(4)) * invs;

tmp1 = q(2) * q(3);
tmp2 = q(4) * q(1);
R(2,1) = 2 * (tmp1 + tmp2) * invs;
R(1,2) = 2 * (tmp1 - tmp2) * invs;

tmp1 = q(2) * q(4);
tmp2 = q(3) * q(1);
R(3,1) = 2 * (tmp1 - tmp2) * invs;
R(1,3) = 2 * (tmp1 + tmp2) * invs;

tmp1 = q(3) * q(4);
tmp2 = q(2) * q(1);
R(3,2) = 2 * (tmp1 + tmp2) * invs;
R(2,3) = 2 * (tmp1 - tmp2) * invs;

end
