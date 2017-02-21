function [ T_ w_ X_ info ] = run_sba( K, T, w, X, x, vis, method, projac_lib )
% RUN_SBA revokes SBA library's bundle adjustment function.
%
% Input:
%       K - camera calibration parameters (4xm)
%       T - camera translation            (3xm)
%       w - camera rotation               (3xm)
%       X - 3d location of points         (4xn)
%       x - image points                  (3xnxm)
%       vis - visibility map              (nxm)
%       method - method of sba optimization ('mot', 'str', or 'motstr')
%       projac_lib - path of projac.dll file
%                    (e.g. '/library/sba-1.4/matlab/projac.dll')
%
% Output:
%       T_ - refined translation
%       w_ - refined rotation
%       X_ - refined 3d location of points
%       info - information from sba()
%
% NOTE: This function is provided for test purpose of SBA library, which
% can be acquired from http://www.ics.forth.gr/~lourakis/sba/

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

% help
if nargin ~= 8
    help run_sba.m;
    return;
end

% get information
m = size(w, 2);
n = size(X, 2);

% exponential map to quaternion
q = zeros(3, m);
for i = 1:m
    qi = iquat( vl_rodr( w(:,i) ) );
    mag = norm(qi);
    if qi(1) < 0
        mag = -mag;
    end
    q(:,i) = qi(2:4) / mag;
end

% prepare sba parameters
p0 = [ reshape([q;T], 1, 6*m) reshape(X(1:3,:), 1, 3*n) ];
cal = [ K(1,1) 0 K(3,1) K(2,1) K(4,1) ];
pts2d = [];
for i=1:n
    for j=1:m
        if vis(i,j)
            pts2d = [ pts2d x(1:2,i,j)' ];
        end
    end
end

% run sba function
% [ret p info] = sba(n, m, 0, vis, p0, 6, 3, pts2d, 2, 'projRTS', 'jacprojRTS', 200, 1, [], method, cal);
imgproj = ['imgproj_motstr@', projac_lib];
imgprojac = ['imgprojac_motstr@', projac_lib];
[ret p info] = sba(n, m, 2, vis, p0, 6, 3, pts2d, 2, imgproj, imgprojac, 200, 1, [], method, cal);

% get refined parameters
cams = reshape(p(1:6*m), 6, m);
q_ = cams(1:3, :);
T_ = cams(4:6, :);
X_ = [ reshape(p(6*m+1:6*m+3*n), 3, n); ones(1,n) ];

% quaternion to exponential map
w_ = zeros(3, m);
for i = 1:m
    qw = sqrt( 1 - sum(q_(:,i).^2) );
    w_(:,i) = vl_irodr( quat( [ qw ; q_(:,i) ] ) );
end

end
