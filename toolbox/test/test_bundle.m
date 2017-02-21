% TEST_BUNDLE tests multiple number of tests on bundle adjustment.

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

num_test = 100;
init_error = zeros(1, num_test);
vlg_error = zeros(1, num_test);
sba_error = zeros(1, num_test);
for i = 1:num_test
    close all;
    disp(['test: ', num2str(i)]);

    % run a single run of bundle adjustment
    demo_bundle_euclid;

    % perform sba library's bundle adjustment
    disp('performing sba ... ');
    [T_sba w_sba X_sba info] = run_sba( K, T1, w1, X1, x, vis, 'motstr', '/library/sba-1.4/matlab/projac.dll' );
    [T_sba w_sba X_sba] = align_scene( T_sba, w_sba, X_sba );
    % display reconstruction
    figure('Name', 'SBA refined scene and motion'), box on, view(3);
    display_points( X, 'pixelsize', 5 );
    display_cameras( T_sba, w_sba, K, 'frustum_scale', 0.1 );
    title('SBA refined scene and motion');
    xlim('auto'), ylim('auto'), zlim('auto');
    % display refinement evaluation
    figure('Name', 'SBA Reconstruction Evaluation');
    subplot(2,3,1), plot(w'),          xlim([1 num_frame]), title('Ground-truth Rotation');
    subplot(2,3,2), plot(w_sba'),      xlim([1 num_frame]), title('SBA Reconstructed Rotation');
    subplot(2,3,3), plot(w' - w_sba'), xlim([1 num_frame]), title('SBA Rotation Error');
    subplot(2,3,4), plot(T'),          xlim([1 num_frame]), title('Ground-truth Translation');
    subplot(2,3,5), plot(T_sba'),      xlim([1 num_frame]), title('SBA Reconstructed Translation');
    subplot(2,3,6), plot(T' - T_sba'), xlim([1 num_frame]), title('SBA Translation Error');
    subplot(2,3,1); yy = ylim; subplot(2,3,2); ylim(yy); drawnow;
    subplot(2,3,4); yy = ylim; subplot(2,3,5); ylim(yy);

    % record the initial and refined error
    init_error(i) = error(1);
    vlg_error(i) = error(end);
    sba_error(i) = info(2)/sum(sum(vis));
end

% plot the histogram of errors
figure('Name', 'Error refinements');
subplot(1,3,1), hist(log10(init_error), 20); title('init error');
subplot(1,3,2), hist(log10(vlg_error), 20); title('vlg');
subplot(1,3,3), hist(log10(sba_error), 20); title('sba');

disp(['mean initial error = ', num2str(mean(init_error))]);
disp(['mean vlg error = ', num2str(mean(vlg_error))]);
disp(['mean sba error = ', num2str(mean(sba_error))]);

