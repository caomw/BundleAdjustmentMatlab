% TEST_MVIEW tests Euclidean 3D reconstruction using multi-view algorithm.
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
number_camera       = 10;   % number of cameras for synthetic data
number_point        = 200;  % number of points for synthetic data
depth               = 100;  % depths of points for synthetic data
pixel_noise         = 0.5;  % pixel noise level for synthetic data

% generate scene and motion, align the reference frame
[K T Omega X x vis] = generate_scene_and_motion( ...
                      number_camera, number_point/2, number_point, depth );

% take points that are visible by all frames for multi-view reconstruction
all_vis_index   = vis(:,1) & vis(:,2);
good_view_index = sum(vis(all_vis_index,:)) > 30;
% select good views from the generated scene
vis = vis(all_vis_index,good_view_index);
X = X(:,all_vis_index);
x = x(:,all_vis_index,good_view_index);
K = K(:,good_view_index);
T = T(:,good_view_index);
Omega = Omega(:,good_view_index);

% align the selected scene
[T Omega X] = align_scene( T, Omega, X );
m = size(x,3);
n = size(x,2);

% add noise to image points
x(1:2,:) = x(1:2,:) + randn(size(x(1:2,:))) * pixel_noise;

% prepare feature tracks
featx(:,:) = x(1,:,:); featx = featx .* vis;
featy(:,:) = x(2,:,:); featy = featy .* vis;

%
% call multi-view reconstruction
%

% % CASE 1 : UNCALIBRATED, USING AN IDENTITY AS A GUESSED CALIBRAITON
% % -----------------------------------------------------------------
% [K_ T_ Omega_ X_ vis_ inlier_] = VLmvg( featx, featy );

% CASE 2 : UNCALIBRATED, BUT PROVIDING A GUESSED CALIBRATION MATRIX
% -----------------------------------------------------------------
% for demo, prepare guessed calibration by adding noise to ground-truth
K1 = K(:,1) + randn(size(K(:,1))) .* 1e-0;
[K_ T_ Omega_ X_ vis_ inlier_] = VLmvg( featx, featy, 'K', calibration_matrix(K1) );
 
% % CASE 3 : CALIBRATED, COORDINATES ARE CALIBRATED ( K = eye(3) )
% % -----------------------------------------------------------------
% [K_ T_ Omega_ X_ vis_ inlier_] = VLmvg( featx, featy, 'calibrated' );

% % CASE 4 : CALIBRATED, PROVIDING A CALIBRATION MATRIX
% % -----------------------------------------------------------------
% [K_ T_ Omega_ X_ vis_ inlier_] = VLmvg( featx, featy, 'calibrated', 'K', calibration_matrix(K(:,1)) );

% % CASE 5 : INITIAL RECONSTRUCTION IS AVAILABLE
% % --------------------------------------------
% % for demo, prepare initial reconstruction by adding noise to ground-truth
% Omega1 = Omega + randn(size(Omega)) * 1e-3;
% T1 = T + randn(size(T)) * 1e-4;
% X1 = X_mv; X1(1:3,:) = X1(1:3,:) + randn(size(X1(1:3,:))) * 1e-3;
% [K_ T_ Omega_ X_ vis_ inlier_] = VLmvg( featx, featy, ...
%                         'calibrated', 'K', calibration_matrix(K(:,1)), ...
%                         'T', T1, 'Omega', Omega1, 'X', X1 );

if sum(inlier_) == 0
    disp('Aborted.');
    return;
end

%
% finishing the bundle adjustment,
%

% compute refined reprojection error
[err_ error_] = error_reproj( x, K_, T_, Omega_, X_, 'visibility', vis_ );
disp(['refined reprojection error = ', num2str(err_)]);
disp([num2str(sum(sign(X_(3,:))>0)), ' positive depths']);
disp([num2str(sum(sign(X_(3,:))<0)), ' negative depths']);
figure('Name', 'Refined reprojection error');
surf(error_);
title('refined reprojection error');
ylabel('features'); xlabel('frames'); zlabel('error');
ylim([1 n]); xlim([1 m]);

% display reconstruction result
figure('Name', 'Euclidean Reconstruction'); box on; view(3);
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
