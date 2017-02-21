function K = calibration_matrix( Kparam )
% CALIBRATION_MATRIX returns a 3x3 calibration matrix from parameters.

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

K = [ Kparam(1) 0         Kparam(3) ; ...
      0         Kparam(2) Kparam(4) ; ...
      0         0         1         ] ;

end