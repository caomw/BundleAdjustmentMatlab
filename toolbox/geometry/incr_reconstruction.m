function [ K T Omega X goodfeat status ] = incr_reconstruction( featx, featy, varargin )
% INCR_RECONSTRUCTION performs incremental reconstruction.
%
% Input:
%        featx - image points x-coordinates (nxm)
%        featy - image points y-coordinates (nxm)
%
% Options:
%        'visibility', vis          - visibility map (nxm)
%        'calibrated'               - indicate that camera is calibrated
%        'K', K (3x3) or (4xm)      - provide a calibration matrix either
%                                     calibration or guessed
%        'initialized', T, Omega, X - provide initial reconstruction for
%                                     the bundle adjustment
%        'status', status           - provide the status of initialized
%                                     frames (1xm logical)
%        'X', X (4xn)               - provide initial points reconstruction
%
% Output:
%        K     - calibration parameters     (4xm)
%        T     - translation                (3xm)
%        Omega - rotation                   (3xm)
%        X     - 3D points                  (4xn) in homogeneous form
%        goodfeat - status of image points  (nxm)
%        status   - status of frames        (1xm) logical

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

% help
num_required_parameters = 2;
if nargin < num_required_parameters
    help incr_reconstruction.m
    return;
end

% get information
m = size(featx, 2);
n = size(featx, 1);

% default values of reconstruction
Omega = zeros(3, m);
T     = zeros(3, m);
X     = zeros(4, n);
status = false(1,m);

% initialize default parameters
vis         = (featx~=0 | featy~=0);
calibrated  = false;
K           = eye(3);
initialized = false;
X_given     = false;

% parse parameters
if nargin > num_required_parameters
    iVarargin = 1;
    while iVarargin <= nargin - num_required_parameters
        switch lower(varargin{iVarargin})
            case 'visibility'
                vis = varargin{iVarargin+1};
                iVarargin = iVarargin + 1;
            case 'calibrated'
                calibrated = true;
            case 'k'
                K = varargin{iVarargin+1};
                iVarargin = iVarargin + 1;
            case 'initialized'
                initialized = true;
                T = varargin{iVarargin+1};
                Omega = varargin{iVarargin+2};
                X = varargin{iVarargin+3};
                iVarargin = iVarargin + 3;
            case 'status'
                status = varargin{iVarargin+1};
                iVarargin = iVarargin + 1;
            case 'x'
                X_given = true;
                X = varargin{iVarargin+1};
                iVarargin = iVarargin + 1;
        end
        iVarargin = iVarargin + 1;
    end
end

% check if we have one K for all images,
% and make calibrations in a 4xm parameter form
if size(K,1) == 3
    % single K for all images
    fix_calibration = true;
    K = repmat( calibration_parameter(K), 1, m );
else
    % different K parameters for each image
    if calibrated
        fix_calibration = true;
    else
        fix_calibration = false;
    end
end

% build measurement matrix x from featx and featy
x = ones(3, n, m);
x(1,:,:) = featx;
x(2,:,:) = featy;

%
% initial reconstruction
%
if initialized == false && X_given == false
    % guess K
    Kguess = calibration_matrix( K(:,1) );

    % choose best pair of images (using eight-point algorithm)
    disp('choosing best two view for initial reconstruction');
    best_i2 = 2;
    best_shared = 0;
    for i2 = 2:m
        % take points that are visible by both frame 1 and i2
        vis_shared = vis(:,1) == 1 & vis(:,i2) == 1;
        [vis_shared_ shared_index] = sort(vis_shared, 'descend');
        n_shared = sum(vis_shared);
        shared_index = shared_index(1:n_shared);

        % take maximum number of shared points
        if n_shared >= min(best_shared, 0.5*sum(vis(:,1)))
            if calibrated
                [Xp lambda] = two_view( x(:,shared_index,1), x(:,shared_index,i2), ...
                    'K', Kguess, 'calibrated' );
            else
                [Xp lambda] = two_view( x(:,shared_index,1), x(:,shared_index,i2), ...
                    'K', Kguess );
            end

            if sum(lambda<=0) == 0
                best_i2 = i2;
                best_shared = n_shared;
            end
        end
    end
    i2 = best_i2;
    disp(['two-view reconstruction between frame 1 and ', num2str(i2)]);

    
    % take points that are visible by both frame 1 and i2
    vis_shared = vis(:,1) == 1 & vis(:,i2) == 1;
    [vis_shared_ shared_index] = sort(vis_shared, 'descend');
    n_shared = sum(vis_shared);
    shared_index = shared_index(1:n_shared);

    % two-view (using eight-point algorithm)
    if calibrated
        [Xp lambda R2 T2] = two_view( ...
            x(:,shared_index,1), x(:,shared_index,i2), ...
            'K', Kguess, 'calibrated' );
    else
        [Xp lambda R2 T2] = two_view( ...
            x(:,shared_index,1), x(:,shared_index,i2), ...
            'K', Kguess );
    end
    
    % in case we have still negative depths, try perturbate focal length
    % until we have all positive depths
    iter_positive_depth = 0;
    while iter_positive_depth < 20 && sum(lambda<=0) > 0
        Kguess_ = Kguess;
        Kguess_(1,1) = Kguess_(1,1)*(1.0+0.01*randn(1));
        Kguess_(2,2) = Kguess_(1,1);
        if calibrated
            [Xp lambda R2 T2] = two_view( ...
                x(:,shared_index,1), x(:,shared_index,i2), ...
                'K', Kguess_, 'calibrated' );
        else
            [Xp lambda R2 T2] = two_view( ...
                x(:,shared_index,1), x(:,shared_index,i2), ...
                'K', Kguess_ );
        end
        iter_positive_depth = iter_positive_depth + 1;
    end

    % if we still have negative depth points, just remove them as outliers
    fprintf('%d negative depth points are removed\n', sum(lambda<0));
    vis(lambda<0,:) = 0;

%     % if we still have negative depth, stop
%     if sum(lambda<=0) > 0
%         K = 0; Omega = 0; T = 0; X = 0; goodfeat = 0;
%         return;
%     end

    % now we have initial reconstruction:
    % Xp, R1 = eye(3), T1 = 0, and R2 and T2
    Omega = zeros(3, m);    Omega(:,i2) = vl_irodr(R2);
    T     = zeros(3, m);    T    (:,i2) = T2;
    X     = zeros(4, n);    X    (:,shared_index)   = Xp;
    status(1)  = true;
    status(i2) = true;
    
    % bundle adjustment with reconstructed points
    X_ba = X(:,shared_index);
    x_ba = x(:,shared_index, [1 i2]);
    [K_ T_ Omega_ X_ba_] = bundle_euclid( ...
        K(:,[1 i2]), T(:,[1 i2]), Omega(:,[1 i2]), X_ba, x_ba, ...
        'fix_calibration', 'verbose' );
    % align the scene and motion to the reference frame and its scale
    [T_ Omega_ X_ba_] = align_scene( T_, Omega_, X_ba_ );
    T    (:,i2)           = T_    (:,2);
    Omega(:,i2)           = Omega_(:,2);
    X    (:,shared_index) = X_ba_;

end % end of initial reconstruction

% incrementally add frames to the reconstruction
% for k = 1:m
%     [dummy j_] = sort( sum(vis(X(4,:) == 1,:), 1), 'descend');
%     j2 = 1;
%     while j2 <= m && status(j_(j2))
%         j2 = j2 + 1;
%     end
%     if j2 > m, continue; end;
%     j = j_(j2);

for j = 1:m
    % skip already added frames
    if status(j), continue; end
    
    disp(['adding image ', num2str(j)]);
    % new camera pose estimation
    vis_index = X(4,:) == 1 & vis(:,j)';
    [K_j T_j Omega_j] = estimate_camera( x(:,vis_index,j), X(:,vis_index), 'K', K(:,j) );
    if size(T_j,1) == 3 && size(Omega_j,1) == 3
        status(j) = true;
    else
        status(j) = false;
        continue;
    end
    
    % add the camera
    T    (:,j) = T_j;
    Omega(:,j) = Omega_j;
    K    (:,j) = K_j;

    % remove outliers
%     [X,vis(:,1:j)] = remove_outlier( K(:,1:j), T(:,1:j), Omega(:,1:j), X, x(:,:,1:j), vis(:,1:j) );
    [X,vis(:,status)] = remove_outlier( K(:,status), T(:,status), Omega(:,status), X, x(:,:,status), vis(:,status) );
    
    % bundle adjustment with reconstructed points
    X3d_index = (X(4,:) == 1);
    X_ba = X(:,X3d_index);
%     x_ba = x(:,X3d_index, 1:j);
%     if fix_calibration
%         [K_ T_ Omega_ X_ba_] = bundle_euclid( ...
%             K(:,1:j), T(:,1:j), Omega(:,1:j), X_ba, x_ba, ...
%             'fix_calibration', 'visibility', vis(X3d_index, 1:j), 'verbose' );
%     else
%         [K_ T_ Omega_ X_ba_] = bundle_euclid( ...
%             K(:,1:j), T(:,1:j), Omega(:,1:j), X_ba, x_ba, ...
%             'visibility', vis(X3d_index, 1:j), 'verbose' );
%     end
    x_ba = x(:,X3d_index, status);
    if fix_calibration
        [K_ T_ Omega_ X_ba_] = bundle_euclid( ...
            K(:,status), T(:,status), Omega(:,status), X_ba, x_ba, ...
            'fix_calibration', 'visibility', vis(X3d_index, status), 'verbose' );
    else
        [K_ T_ Omega_ X_ba_] = bundle_euclid( ...
            K(:,status), T(:,status), Omega(:,status), X_ba, x_ba, ...
            'visibility', vis(X3d_index, status), 'verbose' );
    end
    % align the scene and motion to the reference frame and its scale
    [T_ Omega_ X_ba_] = align_scene( T_, Omega_, X_ba_ );
    % apply refined reconstruction
%     K    (:,1:j) = K_;
%     T    (:,1:j) = T_;
%     Omega(:,1:j) = Omega_;
    K    (:,status) = K_;
    T    (:,status) = T_;
    Omega(:,status) = Omega_;
    X(:,X3d_index) = X_ba_;

    % new 3d point triangulation
    num_new_point = 0;
    for i = 1:n
%         if X(4,i) == 0 && sum(vis(i,1:j)) >= 2
%             vis_i_index = logical(vis(i,1:j));
%             X(:,i) = triangulation( ...
%                         K    (:,vis_i_index), ...
%                         T    (:,vis_i_index), ...
%                         Omega(:,vis_i_index), ...
%                         x  (:,i,vis_i_index) );
        if X(4,i) == 0 && sum(vis(i,status)) >= 2
            vis_i_index = logical(vis(i,:)) & status;
            X(:,i) = triangulation( ...
                        K    (:,vis_i_index), ...
                        T    (:,vis_i_index), ...
                        Omega(:,vis_i_index), ...
                        x  (:,i,vis_i_index) );
            if X(4,i) == 1
                num_new_point = num_new_point + 1;
            end
        end
    end
    disp([num2str(num_new_point), ' new 3d points added.']);

    % remove outliers
%     [X,vis(:,1:j)] = remove_outlier( K(:,1:j), T(:,1:j), Omega(:,1:j), X, x(:,:,1:j), vis(:,1:j) );
    [X,vis(:,status)] = remove_outlier( K(:,status), T(:,status), Omega(:,status), X, x(:,:,status), vis(:,status) );
    
    % bundle adjustment with reconstructed points
    X3d_index = (X(4,:) == 1);
    X_ba = X(:,X3d_index);
%     x_ba = x(:,X3d_index, 1:j);
%     if fix_calibration
%         [K_ T_ Omega_ X_ba_] = bundle_euclid( ...
%             K(:,1:j), T(:,1:j), Omega(:,1:j), X_ba, x_ba, ...
%             'fix_calibration', 'visibility', vis(X3d_index, 1:j), 'verbose' );
%     else
%         [K_ T_ Omega_ X_ba_] = bundle_euclid( ...
%             K(:,1:j), T(:,1:j), Omega(:,1:j), X_ba, x_ba, ...
%             'visibility', vis(X3d_index, 1:j), 'verbose' );
%     end
    x_ba = x(:,X3d_index, status);
    if fix_calibration
        [K_ T_ Omega_ X_ba_] = bundle_euclid( ...
            K(:,status), T(:,status), Omega(:,status), X_ba, x_ba, ...
            'fix_calibration', 'visibility', vis(X3d_index, status), 'verbose' );
    else
        [K_ T_ Omega_ X_ba_] = bundle_euclid( ...
            K(:,status), T(:,status), Omega(:,status), X_ba, x_ba, ...
            'visibility', vis(X3d_index, status), 'verbose' );
    end
    % align the scene and motion to the reference frame and its scale
    [T_ Omega_ X_ba_] = align_scene( T_, Omega_, X_ba_ );
    % apply refined reconstruction
%     K    (:,1:j) = K_;
%     T    (:,1:j) = T_;
%     Omega(:,1:j) = Omega_;
    K    (:,status) = K_;
    T    (:,status) = T_;
    Omega(:,status) = Omega_;
    X(:,X3d_index) = X_ba_;

    figure(10), cla; box on; view(3);
    display_points(X, 'pixelsize', 5);
%     display_cameras(T(:,1:j), Omega(:,1:j), K(:,1:j), 'frustum_scale', 0.1);
    display_cameras(T(:,status), Omega(:,status), K(:,status), 'frustum_scale', 0.1);
    drawnow;
end

goodfeat = vis;

end

function [X vis] = remove_outlier( K, T, Omega, X, x, vis )
    % get information
    m = size(K,2);
    n = size(X,2);
    
    % reproject points
%    error = zeros(n,m);
%    negative = zeros(n,m);
    max_error = 0;
    max_error_point = 0;
    max_error_frame = 0;
    negative_count = 0;
    for j = 1:m
        P_j = calibration_matrix(K(:,j)) * [ vl_rodr(Omega(:,j)) T(:,j) ];
        for i = 1:n
            if X(4,i) == 1 && vis(i,j)
                x_reproj = P_j * X(:,i);
                % check negative depth
                if ( x_reproj(3) < 0 || x_reproj(3) > 10 )
                    negative_count = negative_count + 1;
                    vis(i,j) = 0;
                    X(:,i) = 0;
                else
                    % get reprojection error
                    x_reproj = x_reproj / x_reproj(3);
                    error = norm( x(1:2,i,j) - x_reproj(1:2) )^2;
                    % get max error point
                    if error > max_error
                        max_error = error;
                        max_error_point = i;
                        max_error_frame = j;
                    end
                end
            end
        end
    end

    % remove a point with the largest reprojection error
    if max_error_point > 0 && max_error_frame > 0
        vis(max_error_point, max_error_frame) = 0;
        X(:,max_error_point) = 0;
    end
    
    % remove negative depth points
    if negative_count > 0
        fprintf('%d negative points removed.\n', negative_count);
    end
    
    % remove the farthest point from the center
    max_dist = 0;
    max_i = 0;
    for i = 1:n
        if X(4,i) == 1 && max_dist < norm(X(1:3,i))
            max_dist = norm(X(1:3,i));
            max_i = i;
        end
    end
    if max_i > 0
        X(:,max_i) = 0;
    end
end
% 
% function [X vis] = remove_outlier( K, T, Omega, X, x, vis )
%     % get information
%     m = size(K,2);
%     n = size(X,2);
%     
%     % reproject points
%     error = zeros(n,m);
%     negative = zeros(n,m);
%     for j = 1:m
%         P_j = calibration_matrix(K(:,j)) * [ vl_rodr(Omega(:,j)) T(:,j) ];
%         for i = 1:n
%             if X(4,i) == 1 && vis(i,j)
%                 x_reproj = P_j * X(:,i);
%                 % check negative depth
%                 negative(i,j) = x_reproj(3) < 0;
%                 % get reprojection error
%                 x_reproj = x_reproj / x_reproj(3);
%                 error(i, j) = norm( x(1:2,i,j) - x_reproj(1:2) )^2;
%             end
%         end
%     end
% 
%     % remove a point with the largest reprojection error
%     [error_row max_row] = max(error);
%     [error_col max_col] = max(error_row);
%     vis(max_row(max_col),max_col) = 0;
%     X(:,max_row(max_col)) = 0;
%     
% %     % remove large error points
% %     large_error_points = (error > 100);
% %     if sum(sum(large_error_points)) > 0
% %         fprintf('%d points with large error removed.\n', sum(large_error_points));
% %         vis(large_error_points) = 0;
% %     end
%     
%     % remove negative depth points
%     negative_points = (sum(negative,2)>1);
%     if sum(negative_points) > 0
%         fprintf('%d negative points removed.\n', sum(negative_points));
%         X(:,negative_points) = 0;
%     end
%     
%     % remove the farthest point from the center
%     max_dist = 0;
%     max_i = 0;
%     for i = 1:n
%         if max_dist < norm(X(1:3,i))
%             max_dist = norm(X(1:3,i));
%             max_i = i;
%         end
%     end
%     X(:,max_i) = 0;
% 
% end
% 
