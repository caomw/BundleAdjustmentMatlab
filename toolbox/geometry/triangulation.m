function [ X ] = triangulation( K, T, Omega, x )
% TRIANGULATION computes a 3D location of a point from more than two views.
%
% Input:
%       K     - parameters of calibration matrices         (4xm)
%       T     - camera translations                        (3xm)
%       Omega - camera rotations in exponential map        (3xm)
%       x     - image points in homogeneous representation (3xm)
%
% Output:
%       X     - triangulated 3d location in homogeneous representatin (4x1)
%               , returns X(4) = 0 if any view has negative depth

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

% get information
m = size(K, 2);

% KK = calibration_matrix(K(:,1));
% x = inv(KK) * x;

% build a linear system
A = zeros(2*m, 4);
for j=1:m
    P_j = calibration_matrix(K(:,j)) * [ vl_rodr(Omega(:,j)) T(:,j) ];
%     P_j = [ vl_rodr(Omega(:,j)) T(:,j) ];
    A(2*j-1,:) = x(1,j) * P_j(3,:) - P_j(1,:);
    A(2*j  ,:) = x(2,j) * P_j(3,:) - P_j(2,:);
end

% solve the linear system
[U D V] = svd(A);
X = V(:,4)/V(4,4);

% check positive depth for all views and reprojection error
error = zeros(1,m);
for j=1:m
    P_j = calibration_matrix(K(:,j)) * [ vl_rodr(Omega(:,j)) T(:,j) ];
%     P_j = [ vl_rodr(Omega(:,j)) T(:,j) ];
    x_reproj = P_j * X;
    if x_reproj(3) < 0
        % if any view has negative depth, don't take the triangulation.
        X = zeros(4,1);
        return;
    end
    error(j) = norm( x(1:2,j) - x_reproj(1:2)/x_reproj(3) );
end
% if mean(error) > 1
%     X = zeros(4,1);
% end

end
