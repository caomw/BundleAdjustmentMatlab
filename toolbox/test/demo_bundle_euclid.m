% DEMO_BUNDLE_EUCLID demonstrates the bundle adjustment.

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

% setup parameters
num_frame = 60;
min_feature = 50;
max_feature = 100;
depth = 100;

% generate scene and motion
[K T w X x vis] = generate_scene_and_motion( num_frame, min_feature, max_feature, depth );
[T w X] = align_scene( T, w, X );

figure('Name', 'Ground-truth Scene and motion'); box on; view(3);
display_points( X, 'pixelsize', 5 );
display_cameras( T, w, K, 'frustum_scale', 0.1 );
title('Ground-truth Scene and motion (Synthesized)');
xlim('auto'), ylim('auto'), zlim('auto');

% add noise
% K1 = K + randn(size(K)) * 1e-0;
w1 = w + randn(size(w)) * 1e-3;
T1 = T + randn(size(T)) * 1e-4;
X1 = X; X1(1:3,:) = X1(1:3,:) + randn(size(X1(1:3,:))) * 1e-3;
% x1 = x; x1(1:2,:) = x1(1:2,:) + randn(size(x1(1:2,:))) * 1e-1;

figure('Name', 'Corrupted Scene and motion'); box on; view(3);
display_points( X1, 'pixelsize', 5 );
display_cameras( T1, w1, K, 'frustum_scale', 0.1 );
title('Corrupted Scene and motion');
xlim('auto'), ylim('auto'), zlim('auto');

% perform bundle adjustment
disp('performing bundle adjustment ... ');
tic
[K_ Te_ w_ Xe_ error] = bundle_euclid_nomex( K, T1, w1, X1, x, ...
    'visibility', vis, ...
    'fix_calibration', ...
    'verbose' );
toc

% align to the reference frame and its scale
[Te_ w_ Xe_] = align_scene( Te_, w_, Xe_ );

% plot refined scene and motion
figure('Name', 'VLG refined scene and motion'); box on; view(3);
display_points( Xe_, 'pixelsize', 5 );
display_cameras( Te_, w_, K_, 'frustum_scale', 0.1 );
title('VLG refined scene and motion');
xlim('auto'), ylim('auto'), zlim('auto');

% display refinement evaluation
figure('Name', 'VLG Reconstruction Evaluation');
subplot(2,3,1), plot(w'),        xlim([1 num_frame]), title('Ground-truth Rotation');
subplot(2,3,2), plot(w_'),       xlim([1 num_frame]), title('VLG Reconstructed Rotation');
subplot(2,3,3), plot(w' - w_'),  xlim([1 num_frame]), title('VLG Rotation Error');
subplot(2,3,4), plot(T'),        xlim([1 num_frame]), title('Ground-truth Translation');
subplot(2,3,5), plot(Te_'),      xlim([1 num_frame]), title('VLG Reconstructed Translation');
subplot(2,3,6), plot(T' - Te_'), xlim([1 num_frame]), title('VLG Translation Error');
subplot(2,3,1); yy = ylim; subplot(2,3,2); ylim(yy); drawnow;
subplot(2,3,4); yy = ylim; subplot(2,3,5); ylim(yy);
