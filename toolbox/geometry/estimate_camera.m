function [ K T Omega inlier ] = estimate_camera( x, X, varargin )
% ESTIMATE_CAMERA estimates the camera parameters.
%
% This performs the direct linear transform (DLT) given the set of image
% points x with respect to the corresponding 3D points X. In addition,
% RANSAC alrogithm is used to robustly estimate the camera, followed by the
% non-linear refinement step.
%
% x ~ K [ R | T ] X
%
% input:
%        x = image points (3xn in image plane)
%        X = 3d points (4xn)
%
% Options:
%        'K', K - calibration parameters (4x1) or matrix (3x3)
%
% output:
%        K      - intrinsic camera parameters              (4x1)
%        T      - extrinsic camera parameter (translation) (3x1)
%        Omega  - extrinsic camera parameter (rotation)    (3x1)
%        inlier - inliers among image points

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.
%

% initialize default parameters
calibrated = false;
K = eye(3);

% parse optional parameters
if nargin > 2
    iVarargin = 1;
    while iVarargin <= nargin - 2
        switch lower(varargin{iVarargin})
            case 'k'
                calibrated = true;
                K = varargin{iVarargin+1};
                iVarargin = iVarargin + 1;
        end
        iVarargin = iVarargin + 1;
    end
end

% make K in a matrix form
if size(K,1) == 4
    K = calibration_matrix( K );
end

% get information
n = size(x, 2);
if n < 6
    K      = [];
    T      = [];
    Omega  = [];
    inlier = [];
    disp(['not enough points (', num2str(n), ')']);
    return;
end

% calibrate points if possible
x0 = x;
if calibrated
    x = inv(K) * x0;
end

% normalize x
x_mean = mean(x, 2);
m_x = x_mean(1);
m_y = x_mean(2);
s_x = std(x(1,:));
s_y = std(x(2,:));
Hx = [ 1/s_x 0 -m_x/s_x ; 0 1/s_y -m_y/s_y ; 0 0 1 ];
x_ = Hx * x;

% normalize X
X_mean = mean(X, 2);
m_X = X_mean(1);
m_Y = X_mean(2);
m_Z = X_mean(3);
s_X = std(X(1,:));
s_Y = std(X(2,:));
s_Z = std(X(3,:));
HX = [ 1/s_X 0 0 -m_X/s_X ; ...
       0 1/s_Y 0 -m_Y/s_Y ; ...
       0 0 1/s_Z -m_Z/s_Z ; ...
       0 0 0 1];
X_ = HX * X;

% ransac
ransac_iteration = 100;
num_sample = 6;
dist_threshold = 1e-1;
score_best = 0;
error_best = inf;
for iter = 1:ransac_iteration
    % random sampling
    x_sample = zeros(3,num_sample);
    X_sample = zeros(4,num_sample);
    sample = zeros(n);
    j = 0;
    while j < num_sample
        k = random('unid', n);
        while sample(k) == 1
            k = mod(k, n)+1;
        end
        j = j + 1;
        x_sample(:,j) = x_(:,k);
        X_sample(:,j) = X_(:,k);
        sample(k) = 1;
    end
    
    % perform DLT
    A = zeros(2*num_sample, 12);
    for i = 1:num_sample
        A(2*i-1, :) = [ zeros(1,4), -X_sample(:,i)', x_sample(2,i)*X_sample(:,i)' ];
        A(2*i,   :) = [ X_sample(:,i)', zeros(1,4), -x_sample(1,i)*X_sample(:,i)' ];
    end
    [U D V] = svd(A);
    P_ = reshape(V(:,end), 4, 3)';

    % compute the reprojection error
    error = zeros(1,n);
    x_reproj = P_ * X_;
    for j = 1:n
        if x_reproj(3,j) <= 0
            error(j) = inf;
        else
            error(j) = norm(x_(:,j) - x_reproj(:,j)/x_reproj(3,j));
        end
    end
    
    % consensus
    sample_inlier = error <= dist_threshold;
    score = sum(sample_inlier);
    error_sum = sum(error(sample_inlier));

    % take the best scoring sample set
    if ( score > score_best ) || ...
       ( score == score_best && error_sum < error_best )
        error_best = error_sum;
        score_best = score;
        inlier = sample_inlier;
    end
end
disp(['(inlier, outlier) = ', num2str(sum(inlier)), ', ', num2str(n-sum(inlier)), ...
      ' (error=',num2str(error_best/sum(inlier)), ')']);

if sum(inlier) < 6
    K      = [];
    T      = [];
    Omega  = [];
    disp(['not enough inlier points (', num2str(sum(inlier)), ')']);
    return;
end


% normalize x
x_mean = mean(x(:,inlier), 2);
m_x = x_mean(1);
m_y = x_mean(2);
s_x = std(x(1,inlier));
s_y = std(x(2,inlier));
Hx = [ 1/s_x 0 -m_x/s_x ; 0 1/s_y -m_y/s_y ; 0 0 1 ];
x_ = Hx * x;

% normalize X
X_mean = mean(X(:,inlier), 2);
m_X = X_mean(1);
m_Y = X_mean(2);
m_Z = X_mean(3);
s_X = std(X(1,inlier));
s_Y = std(X(2,inlier));
s_Z = std(X(3,inlier));
HX = [ 1/s_X 0 0 -m_X/s_X ; ...
       0 1/s_Y 0 -m_Y/s_Y ; ...
       0 0 1/s_Z -m_Z/s_Z ; ...
       0 0 0 1];
X_ = HX * X;


% perform DLT only with inliers
A = zeros(2*sum(inlier), 12);
k = 0;
for i = 1:n
    if inlier(i)
        k = k + 1;
        A(2*k-1, :) = [ zeros(1,4), -X_(:,i)', x_(2,i)*X_(:,i)' ];
        A(2*k,   :) = [ X_(:,i)', zeros(1,4), -x_(1,i)*X_(:,i)' ];
    end
end
[U D V] = svd(A);
P_ = reshape(V(:,end), 4, 3)';

% denormalize
P_ = inv(Hx) * P_ * HX;
if calibrated
    P_ = K * P_;
end

% factorize K, R and T
% [K1 R] = rq( P_(1:3,1:3) );
% T = K1 \ P_(1:3,4);
% % normalize K1 to a calibration matrix
% K33_factor = K1(3,3);
% K1 = K1 / K33_factor;
% 
% % make the focal length positive
% if K1(1,1) < 0
%     K1 = K1 * diag([-1,1,1]);
%     R = diag([-1,1,1]) * R;
%     T = diag([-1,1,1]) * T;
% end
% if K1(2,2) < 0
%     K1 = K1 * diag([1,-1,1]);
%     R = diag([1,-1,1]) * R;
%     T = diag([1,-1,1]) * T;
% end
% % make the rotation matrix's determinant = 1
% if det(R) < 0
%     R = -R;
%     T = -T;
% end
if det(P_(1:3,1:3)) < 0
    P_ = -P_;
end
[U D V] = svd(P_(1:3,1:3));
R = U*V';
K1_inv = R*inv(P_(1:3,1:3));
T = K1_inv*P_(1:3,4);

% if K is not calibrated, use K1
if calibrated == false
    K = K1;
end

% build result parameters
K = calibration_parameter( K );
Omega = vl_irodr(R);

%
% non-linear refinement
%
if calibrated
    [K T Omega] = bundle_euclid( K, T, Omega, X, x0, ...
        'fix_calibration', 'fix_structure', 'verbose', 'visibility', inlier' );
else
    [K T Omega] = bundle_euclid( K, T, Omega, X, x0, ...
        'fix_structure', 'verbose', 'visibility', inlier' );
end

end
