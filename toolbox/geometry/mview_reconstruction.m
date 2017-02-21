function [ K T Omega X vis ] = mview_reconstruction( featx, featy, varargin )
% MVIEW_RECONSTRUCTION performs multi-view reconstruction.
%
% Input:
%        featx - image points x-coordinates (nxm)
%        featy - image points y-coordinates (nxm)
%
% Options:
%        'calibrated'               - indicate that camera is calibrated
%        'K', K 3x3                 - provide a calibration matrix either
%                                     calibration or guessed
%        'visibility', vis          - visibility map (nxm)
%        'initialized', T, Omega, X - provide initial reconstruction for
%                                     the bundle adjustment
%        'nomex'                    - use matlab-implementation (no mex)
%                                     of the bundle adjustment
%
% Output:
%        K     - calibration parameters     (4xm)
%        T     - translation                (3xm)
%        Omega - rotation                   (3xm)
%        X     - 3D points                  (4xn)   in homogeneous form
%        vis   - updated visibility map     (nxm)

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

% help
num_required_parameters = 2;
if nargin < num_required_parameters
    help mview_reconstruction.m
    return;
end

% get information
m = size(featx, 2);
n = size(featx, 1);

% initialize default parameters
calibrated = false;
K = eye(3);
initialized = false;
nomex = false;
vis = (featx~=0 | featy~=0);

% parse parameters
if nargin > num_required_parameters
    iVarargin = 1;
    while iVarargin <= nargin - num_required_parameters
        switch lower(varargin{iVarargin})
            case 'calibrated'
                calibrated = true;
            case 'k'
                K = varargin{iVarargin+1};
                iVarargin = iVarargin + 1;
            case 'visibility'
                vis = varargin{iVarargin+1};
                iVarargin = iVarargin + 1;
            case 'initialized'
                initialized = true;
                Te = varargin{iVarargin+1};
                w  = varargin{iVarargin+2};
                Xe = varargin{iVarargin+3};
                iVarargin = iVarargin + 3;
            case 'nomex'
                nomex = true;
        end
        iVarargin = iVarargin + 1;
    end
end

% build measurement matrix x from featx and featy
x = ones(3, n, m);
x(1,:,:) = featx;
x(2,:,:) = featy;

if initialized == false
    % guess K
    Kguess = K;

    % choose best two-view (using eight-point algorithm)
    disp('choosing best two view for initial reconstruction');
    best_i2 = 2;
    for i2 = 2:m
        if sum(vis(:,i2)) == n
            if calibrated
                [Xp lambda] = two_view( x(:,:,1), x(:,:,i2), 'K', Kguess, 'calibrated' );
            else
                [Xp lambda] = two_view( x(:,:,1), x(:,:,i2), 'K', Kguess );
            end

            if sum(lambda<=0) == 0
                best_i2 = i2;
            end
        end
    end

    i2 = best_i2;
    disp(['two-view reconstruction between frame 1 and ', num2str(i2)]);
    if calibrated
        [Xp lambda] = two_view( x(:,:,1), x(:,:,i2), 'K', Kguess, 'calibrated' );
    else
        [Xp lambda] = two_view( x(:,:,1), x(:,:,i2), 'K', Kguess );
    end
    
    % in case we have still negative depths, try perturbate focal length
    % until we have all positive depths
    iter_positive_depth = 0;
    while iter_positive_depth < 20 && sum(lambda<=0) > 0
        Kguess_ = Kguess;
        Kguess_(1,1) = Kguess_(1,1)*(1.0+0.01*randn(1));
        Kguess_(2,2) = Kguess_(1,1);
        if calibrated
            [Xp lambda] = two_view( x(:,:,1), x(:,:,i2), 'K', Kguess_, 'calibrated' );
        else
            [Xp lambda] = two_view( x(:,:,1), x(:,:,i2), 'K', Kguess_ );
        end
        iter_positive_depth = iter_positive_depth + 1;
    end

    % if we still have negative depth points, just remove them as outliers
    fprintf('%d negative depth points are removed\n', sum(lambda<0));
    vis(lambda<0,:) = 0;
    
    % multiple view reconstruction
    disp('multi-view reconstruction');
    disp('-------------------------');
    if calibrated
        [Xp Pp x_ vis] = multi_view( x, lambda, 'K', Kguess, 'calibrated', 'visibility', vis, 'visualize' );
    else
        [Xp Pp x_ vis] = multi_view( x, lambda, 'K', Kguess, 'visibility', vis, 'visualize' );
    end
    if sum(Xp(4,:)) < 10
        disp('multi-view reconstruction failed');
        K = 0; Omega = 0; T = 0; X = 0;
        return;
    end
    
    if calibrated == false
        % in case uncalibrated case, perform projective bundle adjustment
        % followed by euclidean upgrade

        % projective bundle adjustment
        disp('projective bundle adjustment');
        disp('----------------------------');
        [Pp_ Xp_] = bundle_projective( Pp, Xp, x_, 'visibility', vis, 'verbose' );

        % euclidean upgrade
        disp('euclidean upgrade');
        disp('-----------------');
        [Xe Re Te K] = upgrade_euclid( Xp_, Pp_ );
        
        % w (exponential map)
        w = zeros(3,m);
        for i=1:m
            w(:,i) = vl_irodr( Re(:,:,i) );
        end
    else
        % in case calibrated case, directly go to euclidean bundle adjustment

        % take the projective reconstruction as euclidean (ideally)
        Xe = Xp;
        Te(:,:) = Pp(:,4,:);
        % w (exponential map)
        w = zeros(3,m);
        for i=1:m
            w(:,i) = vl_irodr( Pp(1:3,1:3,i) );
        end

        % take the initial K as calibrated ones for all frames
        K = zeros(3,3,m);
        for i = 1:m
            K(:,:,i) = eye(3,3);
        end
    end

    % assume fixed camera parameters over all frames
    K = Kguess * K(:,:,1);
    disp('K = ');
    disp(K);
    % reprojection error
    err = error_reproj( x, K, Te, w, Xe, 'visibility', vis );
    disp(['(assume fixed camera parameters) reprojection error = ', num2str(err)]);

end % end of initial reconstruction

% build calibration parameter matrix
Kparam = repmat(calibration_parameter(K), 1, m);

% euclidean bundle adjustment
disp('euclidean bundle adjustment');
disp('---------------------------');
if nomex == false
    [K T Omega X] = bundle_euclid( Kparam, Te, w, Xe, x, 'fix_calibration', 'visibility', vis, 'verbose' );
else
    [K T Omega X] = bundle_euclid_nomex( Kparam, Te, w, Xe, x, 'fix_calibration', 'visibility', vis, 'verbose' );
end

end
