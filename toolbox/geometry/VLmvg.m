function [ K_ T_ Omega_ X_ vis_ inlier_ ] = VLmvg( featx, featy, varargin )
% VLMVG reconstructs structure and motion from feature measurements.
%
% Input:
%        featx - image points x-coordinates (nxm)
%        featy - image points y-coordinates (nxm)
%
% Options:
%        'visibility', vis - visibility map     (nxm)
%        'calibrated'      - indicates that the camera is calibrated
%        'K', K            - camera calibraiton (3x3 or 4xm)
%        'T', T            - camera translation (3xm)
%        'Omega', Omega    - camera rotation    (3xm) in exponential map
%        'X', X            - point location     (4xn) in homogeneous form
%        'ransac'          - performs ransac for rejecting outliers
%        'sampson_threshold', threshold - sampson threshold (default=1)
%
% Output:
%        K_       - calibration parameters  (4xm)
%        T_       - camera translation      (3xm)
%        Omega_   - camera rotation         (3xm)
%        X_       - 3D point location       (4xn) in homogeneous form
%        vis_     - status of image points  (nxm)
%        inlier_  - indicates inlier points (1xn)

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

% help
num_required_parameters = 2;
if nargin < num_required_parameters
    help VLmvg.m
    return;
end

% get information
m = size(featx, 2);
n = size(featx, 1);

% initialize default parameters
vis         = (featx~=0 | featy~=0);
calibrated  = false;
K_given     = false;
T_given     = false;
Omega_given = false;
X_given     = false;
perform_ransac    = false;
sampson_threshold = 1;

% parse parameters
try
    if nargin > num_required_parameters
        iVarargin = 1;
        % for backward compatibility, we specially check the visibility
        if ndims(varargin{1}) == 2
            if size(varargin{1},1) == n && size(varargin{1},2) == m
                vis = varargin{1};
                iVarargin = iVarargin + 1;
            end
        end
        % parse the optional parameters
        while iVarargin <= nargin - num_required_parameters
            switch lower(varargin{iVarargin})
                case 'visibility'
                    vis = varargin{iVarargin+1};
                    iVarargin = iVarargin + 1;
                case 'calibrated'
                    calibrated = true;
                case 'k'
                    K_given = true;
                    K = varargin{iVarargin+1};
                    iVarargin = iVarargin + 1;
                case 't'
                    T_given = true;
                    T = varargin{iVarargin+1};
                    iVarargin = iVarargin + 1;
                case 'omega'
                    Omega_given = true;
                    Omega = varargin{iVarargin+1};
                    iVarargin = iVarargin + 1;
                case 'x'
                    X_given = true;
                    X = varargin{iVarargin+1};
                    iVarargin = iVarargin + 1;
                case 'ransac'
                    perform_ransac = true;
                case 'sampson_threshold'
                    sampson_threshold = varargin{iVarargin+1};
                    iVarargin = iVarargin + 1;
            end
            iVarargin = iVarargin + 1;
        end
    end
catch
    error('Invalid parameters.');
end

% now that we have option parameters given, we choose an appropriate
% algorithm to run the reconstruction pipeline.

if perform_ransac == true
    % in case the feature track is not polished, we reject outliers by
    % applying epipolar constraints between frames.
    vis = remove_outliers( featx, featy, vis, 'threshold', sampson_threshold );
end

if calibrated
    % in case camera is calibrated
    
    if K_given == false
        % in this case, we know that the xfeat and yfeat are calibrated,
        % and K is an identity matrix.
        K = eye(3);
    end
    
    if T_given && Omega_given && X_given
        % in case initial reconstruction is given,
        if size(K,1) == 3, K = repmat(calibration_parameter(K), 1, m); end
        x = ones(3, n, m);
        x(1,:,:) = featx;
        x(2,:,:) = featy;
        % euclidean bundle adjustment
        disp('euclidean bundle adjustment');
        disp('---------------------------');
        [K_ T_ Omega_ X_] = bundle_euclid( K, T, Omega, X, x, 'visibility', vis, 'fix_calibration', 'verbose' );
        vis_ = vis;
    else
        % in case there is no initial reconstruction,
        if sum(vis(:,1)) == n && sum(vis(:,2)) == n
            % in case all features are visible by all frames
            [K_ T_ Omega_ X_ vis_] = mview_reconstruction( featx, featy, ...
                'calibrated', 'K', K , 'visibility', vis );
        else
            % in case some features are not visible by some frames
            [K_ T_ Omega_ X_ vis_] = incr_reconstruction( featx, featy, ...
                'calibrated', 'K', K, 'visibility', vis );
        end
    end
else
    % camera is not calibrated
    
    if K_given == false
        % in case no information on K is given, we guess one from features
        cx = ( max(max(featx)) - min(min(featx)) ) / 2;
        cy = ( max(max(featy)) - min(min(featy)) ) / 2;
        fx = 2 * cx - min(min(featx));
        fy = fx;
        K  = [ fx 0 cx ; 0 fy cy ; 0 0 1 ];
    end

    if T_given && Omega_given && X_given
        % in case initial reconstruction is given,
        if size(K,1) == 3, K = repmat(calibration_parameter(K), 1, m); end
        x = ones(3, n, m);
        x(1,:,:) = featx;
        x(2,:,:) = featy;
        % euclidean bundle adjustment
        disp('euclidean bundle adjustment');
        disp('---------------------------');
        [K_ T_ Omega_ X_] = bundle_euclid( K, T, Omega, X, x, 'visibility', vis, 'verbose' );
        vis_ = vis;
    else
        if sum(vis(:,1)) == n && sum(vis(:,2)) == n
            % in case all features are visible by all frames
            [K_ T_ Omega_ X_ vis_] = mview_reconstruction( featx, featy, ...
                'K', K, 'visibility', vis );
        else
            % in case some features are not visible by some frames
            [K_ T_ Omega_ X_ vis_] = incr_reconstruction( featx, featy, ...
                'K', K, 'visibility', vis );
        end
    end
end

% align the scene
[T_ Omega_ X_] = align_scene( T_, Omega_, X_ );

% inliers are the points whose 3d locations are successfully reconstructed.
inlier_ = X_(4,:) == 1;

end
