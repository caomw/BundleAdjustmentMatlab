function q = iquat( R )
% IQUAT computes a quaternion from a rotation matrix R
%
% Input:
%       R - 3x3 rotation matrix
%
% Output:
%       q - 4x1 quaternion vector

% reference: Conversion Matrix to Quaternion - Martin Baker

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

tr = trace(R) + 1;
if tr > eps
    s = 0.5 / sqrt(tr);
    w = 0.25 / s;
    x = (R(3,2) - R(2,3)) * s;
    y = (R(1,3) - R(3,1)) * s;
    z = (R(2,1) - R(1,2)) * s;
else
    if R(1,1) > R(2,2) && R(1,1) > R(3,3)
        s = 2 * sqrt( 1 + R(1,1) - R(2,2) - R(3,3) );
        w = (R(2,3) - R(3,2)) / s;
        x = 0.25 * s;
        y = (R(1,2) + R(2,1)) / s;
        z = (R(1,3) + R(3,1)) / s;
    elseif R(2,2) > R(3,3)
        s = 2 * sqrt( 1 + R(2,2) - R(1,1) - R(3,3) );
        w = (R(1,3) - R(3,1)) / s;
        x = (R(1,2) + R(2,1)) / s;
        y = 0.25 * s;
        z = (R(2,3) + R(3,2)) / s;
    else
        s = 2 * sqrt( 1 + R(3,3) - R(1,1) - R(2,2) );
        w = (R(1,2) - R(2,1)) / s;
        x = (R(1,3) + R(3,1)) / s;
        y = (R(2,3) + R(3,2)) / s;
        z = 0.25 * s;
    end
end

q = [ w x y z ]';

end
