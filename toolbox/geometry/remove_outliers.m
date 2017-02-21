function [vis_] = remove_outliers( featx, featy, vis, varargin )
% REMOVE_OUTLIERS removes outliers among feature measurements.
%
% Input:
%        featx - image points x-coordinates (nxm)
%        featy - image points y-coordinates (nxm)
%        vis   - visibility map             (nxm)
%
% Options:
%        'iteration', iteration - ransac iteration
%        'threshold', threshold - sampson distance threshold
%
% Output:
%        vis_ - visibility map with inliers (nxm)

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

% help
num_required_parameters = 3;
if nargin < num_required_parameters
    help remove_outliers.m
    return;
end

% initialize default parameters
ransac_iteration    = 100;  % number of iterations for ransac algorithm
sampson_threshold   = 1e-6; % sampson distance threshold

% parse parameters
try
    if nargin > num_required_parameters
        iVarargin = 1;
        while iVarargin <= nargin - num_required_parameters
            switch lower(varargin{iVarargin})
                case 'iteration'
                    ransac_iteration = varargin{iVarargin+1};
                    iVarargin = iVarargin + 1;
                case 'threshold'
                    sampson_threshold = varargin{iVarargin+1};
                    iVarargin = iVarargin + 1;
            end
            iVarargin = iVarargin + 1;
        end
    end
catch
    error('Invalid parameters.');
end

% get information
m = size(featx, 2);

% ransac
disp('ransac with epipolar constraint');
disp('-------------------------------');
disp(['sampson distance threshold = ', num2str(sampson_threshold)]);
for i1 = 1:m-1
    for i2 = i1+1:m
        % prepare feature measurements for each frame
        tracked = vis(:,i1) == 1 & vis(:,i2) == 1;
        n_tracked = sum(tracked);
        if n_tracked > 30
            [tracked_ index] = sort(tracked, 'descend');
            index = index(1:n_tracked);
            x1 = [ featx(index,i1)' ; featy(index,i1)' ; ones(1,n_tracked) ];
            x2 = [ featx(index,i2)' ; featy(index,i2)' ; ones(1,n_tracked) ];

            % apply epipolar constraint
            inlier = ransac_epipolar_constraint( ...
                x1, x2, ransac_iteration, sampson_threshold );

            % update the visibility map
            vis(index,i2) = vis(index,i2) .* inlier';
        end
    end
end

% return updated visibility map
vis_ = vis;

end

