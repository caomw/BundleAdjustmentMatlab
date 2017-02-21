% TEST_INCREMENTAL tests Euclidean reconstruction using incremental bundle adjustment.
%
% Taehee Lee
%
% last update: 06/16/2008
%

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

close all;
clear;

%
% set parameters for test
%
number_camera       = 50;  % number of cameras for synthetic data
number_point        = 200;  % number of points for synthetic data
depth               = 100;  % depths of points for synthetic data
pixel_noise         = 0.5;  % pixel noise level for synthetic data

% generate scene and motion, align the reference frame
[K T Omega X x vis] = generate_scene_and_motion( ...
                      number_camera, number_point/2, number_point, depth );
[T Omega X] = align_scene( T, Omega, X );

m = number_camera;

% add noise to image points
x(1:2,:) = x(1:2,:) + randn(size(x(1:2,:))) * pixel_noise;

% prepare feature tracks
featx(:,:) = x(1,:,:); featx = featx .* vis;
featy(:,:) = x(2,:,:); featy = featy .* vis;

%
% call incremental reconstruction
%

tic
[ K_ T_ Omega_ X_ goodfeat inlier_ ] = VLmvg( featx, featy, 'calibrated', 'K', K );
toc

if sum(inlier_) == 0
    disp('Aborted.');
    return;
end

%
% finishing the bundle adjustment,
%

% get reconstructed number of points
n = size(X_,2);

% compute refined reprojection error
[err_ error_] = error_reproj( x, K_, T_, Omega_, X_, 'visibility', goodfeat );
disp(['refined reprojection error = ', num2str(err_)]);
figure('Name', 'Refined reprojection error');
surf(error_);
title('refined reprojection error');
ylabel('features'); xlabel('frames'); zlabel('error');
ylim([1 n]); xlim([1 m]);

% display reconstruction result
figure('Name', 'Euclidean Reconstruction'), box on, view(3);
display_points( X,  'pixelsize', 5, 'color', [1 0 0] ); % ground-truth
display_points( X_, 'pixelsize', 5, 'color', [0 0 1] ); % reconstruction
% display ground-truth and reconstructed camera trajectory
display_cameras( T,  Omega,  K,  'no_frustum', 'trajectory_color', [1 0 0 ]);
display_cameras( T_, Omega_, K_, 'frustum_scale', 0.1 );
title('Euclidean Reconstruction (Red: Ground-truth, Blue: Reconstruction)');
xlim('auto'), ylim('auto'), zlim('auto');

% display refinement evaluation
figure('Name', 'Reconstruction Evaluation');
subplot(2,3,1), plot(Omega'), xlim([1 m]), title('Ground-truth Rotation');
subplot(2,3,2), plot(Omega_'), xlim([1 m]), title('Reconstructed Rotation');
subplot(2,3,3), plot(Omega' - Omega_'), xlim([1 m]), title('Rotation Error');
subplot(2,3,4), plot(T'), xlim([1 m]), title('Ground-truth Translation');
subplot(2,3,5), plot(T_'), xlim([1 m]), title('Reconstructed Translation');
subplot(2,3,6), plot(T' - T_'), xlim([1 m]), title('Translation Error');
subplot(2,3,1); yy = ylim; subplot(2,3,2); ylim(yy); drawnow;
subplot(2,3,4); yy = ylim; subplot(2,3,5); ylim(yy);
