% DEMO_VLMVG demonstrates an example of 3D reconstruction pipeline.

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

% change the directory to the vlg/data
currentpath = pwd;
cd([vlg_root,'/data']);

% get a list of image files
filelist = dir('UCLA*.jpg');

% get feature track from the image files
[featx,featy,vis,color] = generate_feature_track( filelist );

% run mview reconstruction
% index = sum(vis,2) == size(vis,2); % take features that are seen by all views
index = vis(:,1) & vis(:,2); % take features that are seen by first and second views
good_view_index = sum(vis(index,:)) > 30;
[K T Om X] = VLmvg(featx(index,good_view_index), featy(index,good_view_index) );
% [K T Om X] = VLmvg(featx(index,good_view_index), featy(index,good_view_index), vis(index,good_view_index) );
figure, display_points(X, 'pixelsize', 5, 'color', color(index,:)), display_cameras(T, Om, K, 'frustum_scale', 0.1), box on, view(3);

%
% Instead of mview reconstruction, you may want to run incremental reconstruction.
%
% % run incremental reconstruction
% [K T Om X] = VLmvg(featx, featy, vis );
% figure, display_points(X, 'pixelsize', 5, 'color', color),
% display_cameras(T, Om, K, 'frustum_scale', 0.1), box on, view(3);

% display image pair with correspondences
i1 = 1;
i2 = 2;
display_image_pair(imread(filelist(i1).name), imread(filelist(i2).name), [featx(:,i1) featy(:,i1)]', [featx(:,i2) featy(:,i2)]');

% return to the original path
cd(currentpath);
clear currentpath;
