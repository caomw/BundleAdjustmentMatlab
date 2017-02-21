function [ featx featy vis color ] = generate_feature_track( filelist )
% GENERATE_FEATURE_TRACK generates a feature track from a list of images.
%
% Input:
%        filelist - (mx1) struct array, containing filenames of images with
%                   its field 'name', which is the output of 'dir' command.
%
% Output:
%        featx - image points x-coordinates (nxm)
%        featy - image points y-coordinates (nxm)
%        vis   - visibility map             (nxm)
%        color - color of points in RGB     (nx3)

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

% get information
m = size(filelist, 1);

% detect sift features
f = cell(m,1);  % sift keypoints
d = cell(m,1);  % sift descriptors
c = cell(m,1);  % colors
for i = 1:m
    fprintf('detecting sift keypoints from %s ... ', filelist(i).name);
    I0 = imread(filelist(i).name);
    if size(I0,3) == 3
        I = rgb2gray(I0);
    else
        I = I0;
    end
    [f{i} d{i}] = vl_sift(single(I));
    if size(I0,3) == 3
        c{i} = zeros(size(f{i},2),3);
        for j = 1:size(f{i},2)
            y = min(max(ceil(f{i}(2,j)),1),size(I0,1));
            x = min(max(ceil(f{i}(1,j)),1),size(I0,2));
            c{i}(j,:) = reshape(I0(y,x,:), 1, 3);
        end
    else
        c{i} = zeros(size(f{i},2),3);
        for j = 1:size(f{i},2)
            y = min(max(ceil(f{i}(2,j)),1),size(I,1));
            x = min(max(ceil(f{i}(1,j)),1),size(I,2));
            c{i}(j,:) = repmat(I(y,x), 1, 3);
        end
    end
    fprintf('%d found.\n', size(f{i},2));
end

% construct a feature track from the first frame
N = size(f{1}, 2);
Track(1).ID = 1:N;
Track(1).x = [ f{1}(1:2,:) ; ones(1,N) ];
Desc = d{1};
color = c{1};

% match other frames
for i = 2:m
    fprintf('frame %d:\n', i);
    
    % initialize the Track of frame i
    Track(i).ID = zeros(1,size(f{i},2)); %#ok<AGROW>
    Track(i).x = [ f{i}(1:2,:) ; ones(1,size(f{i},2)) ]; %#ok<AGROW>
    
    % match the features to other frames
    disagree = zeros(1,size(f{i},2));
     for j = 1:i-1 % matching all frames to each other
%      for j = max(1,i-1):i-1 % matching between most current two frames
%     for j = max(1,i-2):i-1 % matching among most current three frames
        % check if we have matched all features
        if sum(Track(i).ID == 0) == 0
            break;
        end
        
        % match between frame i and j
        matches = match_sift_unique( d{i}, Desc(:,Track(j).ID) );
        x1 = [ f{i}((1:2),matches(1,:)) ; ones(1,size(matches,2)) ];
        x2 = Track(j).x(:,matches(2,:));
        inlier = ransac_epipolar_constraint(x1, x2, 100, 1);
        
        % add the matched features
        if sum(inlier) > 30
            old_id = Track(i).ID(matches(1,inlier));
            new_id = Track(j).ID(matches(2,inlier));
            disagree(matches(1,old_id ~= new_id & old_id ~= 0)) = 1;
            Track(i).ID(matches(1,inlier)) = Track(j).ID(matches(2,inlier)); %#ok<AGROW>
        end
    end
    Track(i).ID(logical(disagree)) = 0; %#ok<AGROW>
    
    % add the unmatched features
    new_feature_index = Track(i).ID == 0;
    num_new_features = sum(new_feature_index);
    Track(i).ID(new_feature_index) = N+1:N+num_new_features; %#ok<AGROW>
    Desc(:,N+1:N+num_new_features) = d{i}(:,new_feature_index);
    color(N+1:N+num_new_features,:) = c{i}(new_feature_index,:);
    N = N + num_new_features;
end

% construct output parameters
featx = zeros(N,m);
featy = zeros(N,m);
vis = zeros(N,m);
for i = 1:m
    featx(Track(i).ID,i) = Track(i).x(1,:);
    featy(Track(i).ID,i) = Track(i).x(2,:);
    vis(Track(i).ID,i) = Track(i).x(3,:);
end

% take only the points that are visible by more than two views
good = sum(vis,2) > 2;
vis = logical(vis(good,:));
featx = featx(good,:); featx = featx .* vis;
featy = featy(good,:); featy = featy .* vis;
color = color(good,:) ./ max(max(color));

end
