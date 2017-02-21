function [r q] = rq( A )
% RQ performs RQ-decomposition of A.
%
% Input:
%        A - input matrix
%
% Output:
%        r - upper triangular matrix
%        q - orthogonal matrix
%        such that A = r q

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

[ Et Lt ] = qr( inv(A) );
r = inv(Lt);
q = inv(Et);

end
