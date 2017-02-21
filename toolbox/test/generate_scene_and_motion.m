function [K T w X x vis] = generate_scene_and_motion( m, min_n, max_n, depth )
% GENERATE_SCENE_AND_MOTION generates a synthetic scene and camera motion.
%
% Input:
%       m - number of frames
%       min_n - minimum number of features to track in a frame
%       max_n - maximum number of features to track in a frame
%       depth - mean depth of points to generate at a certain frame
%
% Output:
%       K - camera calibration parameters (4xm)
%       T - camera translation            (3xm)
%       w - camera rotation               (3xm)
%       X - 3d location of points         (4xn)
%       x - image points                  (3xnxm)
%       vis - visibility map              (nxm)

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

% help
if nargin ~= 4
    help generate_scene_and_motion.m;
    return;
end

H = figure('Name', 'Generated Tracks');

% generate image parameters
width  = 500;
height = 500;

% generate calibration parameters
f  = 1*width;       % focal length in pixels (assume square pixels)
cx = width  / 2;    % principal point's x coordinate
cy = height / 2;    % principal point's y cooridnate
K = repmat( [ f f cx cy ]', 1, m );

% generate camera motion
w = zeros(3,m);
T = zeros(3,m);
vel_w = zeros(3,1);
vel_T = zeros(3,1);
for j=2:m
    acc_w = 1e-2*randn(3,1);    % random rotational acceleration model
    acc_T = 1e-0*randn(3,1);    % random translational acceleration model
    vel_w = vel_w + acc_w;
    vel_T = vel_T + acc_T;
    w(:,j) = w(:,j-1) + vel_w;
    T(:,j) = T(:,j-1) + vel_T;
end

% generate measurements
X   = zeros( 4, max_n );        % the scene points
x   = zeros( 3, max_n, m );     % this measurement matrix will grow
vis = zeros(max_n,m); % and also this visibility map will grow
n   = 0;
for j=1:m
    % track previous features
    n_tracked = 0;
    for i=1:n
        if vis(i,j-1) == 1
            % project the 3d point to the current frame
            x_ij = calibration_matrix(K(:,j)) * [ vl_rodr(w(:,j)) T(:,j) ] * X(:,i);
            if x_ij(3) > 0.01*depth && ...
               x_ij(1)/x_ij(3) > 1 && x_ij(1)/x_ij(3) < width && ...
               x_ij(2)/x_ij(3) > 1 && x_ij(2)/x_ij(3) < height
                % if the projection is visible at this frame, take it
                n_tracked = n_tracked + 1;
                x(:,i,j) = x_ij / x_ij(3);
                vis(i,j) = 1;
            else
                x(:,i,j) = 0;
                vis(i,j) = 0;
            end
        end
    end
    % detect previously stored features
    if n_tracked < min_n
        for i=1:n
            if vis(i,j) == 0
                % project the 3d point to the current frame
                x_ij = calibration_matrix(K(:,j)) * [ vl_rodr(w(:,j)) T(:,j) ] * X(:,i);
                if n_tracked < max_n && ...
                   x_ij(3) > 0.01*depth && ...
                   x_ij(1)/x_ij(3) > 1 && x_ij(1)/x_ij(3) < width && ...
                   x_ij(2)/x_ij(3) > 1 && x_ij(2)/x_ij(3) < height
                    % if the projection is visible at this frame, take it
                    n_tracked = n_tracked + 1;
                    x(:,i,j) = x_ij / x_ij(3);
                    vis(i,j) = 1;
                end
            end
        end
    end
    % introduce new features
    if n_tracked < min_n
        n_new = max_n - n_tracked;
        % generate features in the current image plane
        xi = rand(2,n_new) .* repmat([ width height ]', 1, n_new);
        d  = (1 + 0.5 * randn(1,n_new))  * depth;
        % then compute the 3d locations
        K_inv = pinv( calibration_matrix(K(:,j)) );
        for i = 1:n_new
            if d(i) > 0
                Xi = vl_rodr(w(:,j))' * ( d(i) * K_inv * [ xi(:,i) ; 1 ] - T(:,j) );
                % add to scene points
                n = n + 1;
                X(:,n) = [ Xi ; 1 ];
                % add to measurements
                x(:,n,j) = [ xi(:,i) ; 1 ];
                vis(n,j) = 1;
            end
        end
    end
    
    figure(H), imagesc(~vis), colormap(gray);
end

end
