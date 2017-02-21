function Kparam = calibration_parameter( K )
% CALIBRATION_PARAMETER returns calibration parameters (4x1) from a matrix.

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

Kparam = [ K(1,1) K(2,2) K(1,3) K(2,3) ]';

end