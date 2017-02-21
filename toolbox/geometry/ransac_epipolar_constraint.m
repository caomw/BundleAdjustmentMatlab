function [ inlier ] = ransac_epipolar_constraint( x1, x2, iter, thresh )
% RANSAC_EPIPOLAR_CONSTRAINT applies the epipolar constraint using RANSAC.
%
% Input:
%       x1, x2 = 3-by-n image point measurement for two frames
%       iter   = maximum number of iteration
%       thresh = threshold of the sampson distance
%
% Output:
%       inlier = 1-by-n array indicating inliers
%                (e.g. [ 0 1 1 1 1 ] --> the first element is an outlier.)

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

% get information
n = size(x1, 2);

% RANSAC
inlier = ones(1, n);
ransac_iteration = iter;
num_sample = 10;
sampson_threshold = thresh;
score_best = 0;
sampson_best = inf;

% check if we have enough number of samples
if n < 2 * num_sample
    % not enough samples, thus return all as inliers
    disp(['not enough samples (', num2str(n), ')']);
    inlier = ones(1,n);
else
    % enough samples, do the ransac
    
    for i=1:ransac_iteration
        x1_ = zeros(3, num_sample);
        x2_ = zeros(3, num_sample);
        % random sample
        sample = zeros(n);
        j = 0;
        while j < num_sample
            k = random('unid', n);
            if sample(k) == 0
                j = j + 1;
                x1_(:,j) = x1(:,k);
                x2_(:,j) = x2(:,k);
                sample(k) = 1;
            end
        end

        % compute F
        F = fundamental( x1_, x2_ );

        % compute sampson distance
        sampson = sampson_distance( x1, x2, F );

        % consensus
        sample_inlier = sampson < sampson_threshold;
        score = sum(sample_inlier);
        sampson_sum = sample_inlier' * sampson ;

        % take the best scoring sample set
        if ( score > score_best ) || ...
           ( score == score_best && sampson_sum < sampson_best )
            sampson_best = sampson_sum;
            score_best = score;
            inlier = sample_inlier';
        end
    end
    disp(['(inlier, outlier) = ', num2str(sum(inlier)), ', ', num2str(n-sum(inlier))]);
end

end
