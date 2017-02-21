function d = sampson_distance( x1, x2, F )
% SAMPSON_DISTANCE computes the sampson distance between x1 and x2 given F
%
% Input:
%       x1, x2 : image points of two frames (3xn)
%       F      : fundamental matrix between the two frames (3x3)
%
% Output:
%       d - the sampson distance

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

n = size(x1, 2);
e3 = [ 0 ; 0 ; 1 ];
e3hat = vl_hat(e3);

d = zeros(n,1);
for i = 1:n
    d(i) = ( x2(:,i)' * F * x1(:,i) )^2 ...
         / ( norm(e3hat * F * x1(:,i))^2 + ...
             norm(x2(:,i)' * F * e3hat)^2 );
end

end
