function [ Xp Pp x_ vis ] = multi_view( x, depth, varargin )
% MULTI_VIEW performs multi-view projective reconstruction.
%
% Input:
%        x     - image points (3xnxm)
%        depth - initial depths of points from two-view (1xn)
%
% Options:
%        'K', K            - calibration matrix (3x3)
%        'calibrated'      - indicates the camera is calibrated
%        'visibility', vis - visibility map     (nxm)
%        'visualize'       - visualize intermediate reconstruction status
%
% Output:
%        Xp  - 3d locations of points  (4xn)
%        Pp  - projection matrices     (3x4xm)
%        x_  - normalized image points (3xnxm)
%        vis - updated visibility map  (nxm)
%
% NOTE: We assume tvl_hat the points in the first image are never outliers.

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

% help
num_required_parameters = 2;
if nargin < num_required_parameters
    help multi_view.m
    return;
end

% get information
m  = size(x,3);
n  = size(x,2);

% initialize default parameters
calibrated = false;
K = eye(3);
vis = ones(n,m);
visualize = false;

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
            case 'visualize'
                visualize = true;
        end
        iVarargin = iVarargin + 1;
    end
end

% open a status figure
if visualize == true
    fig_status = figure(10);
    screen_size = get(0, 'ScreenSize');
    fig_w = 1600;
    fig_h = 400;
    set(fig_status, 'Position', [ 0 screen_size(4)-fig_h-100 fig_w fig_h ], ...
                    'Name', 'Multi-view reconstruction status');
end

% normalize x by (guessed) K
Kinv = inv(K);
x_ = zeros(3,n,m);
for i = 1:m
    x_(:,:,i) = Kinv * x(:,:,i);
end

% init R and T for reference frame
R(:,:,1) = eye(3);
T(:,1) = zeros(3,1);
Pp(:,:,1) = [ R(:,:,1) T(:,1) ];

% initialize the scale
alpha = 1 ./ depth;
% scale = alpha(1);
scale = mean(alpha(logical(vis(:,1))));

% init X with initial depths
Xp = zeros(4,n);
for i=1:n
    if alpha(i) > 0 && vis(i,1)
        Xp(:,i) = [scale/alpha(i)*x_(:,i,1) ; 1];
    end
end

% check if original projective depths seem to be all positive
% all_positive_depth = sum(sign(alpha)<0) < 0.05*n; % less than 5% negatives
all_positive_depth = 1;

%
% loop until error converges, or until specified iterations
%
old_err = inf;
err = old_err;
iter = 0;
err_increase_iter = 0;

while (iter<2 || abs(old_err-err)>1e-4*old_err) && ...
       iter<20 && err_increase_iter < 3
    %
    % 2. set alpha
    %
    alpha = alpha / scale;

    %
    % 3. assemble Pi and compute Ri and Ti
    %
    R_ = zeros(3,3,m);
    T_ = zeros(3,m);
    R_(:,:,1) = eye(3);
    T_(:,1) = 0;
    for i = 2:m
        % assemble Pi by taking visible points
        Pi = zeros(3*sum(vis(:,i)),12); % Pi = (3n-by-12)
        k = 0;
        for j = 1:n
            if vis(j,i)
                k = k + 1;
                Pi(3*k-2:3*k,:) = [ kron(x_(:,j,1)',vl_hat(x_(:,j,i))) alpha(j)*vl_hat(x_(:,j,i)) ];
            end
        end

        % compute Ri and Ti
        [U D V] = svd(Pi);
        R_(:,:,i) = reshape( V(1:9, 12), 3, 3 );
        T_(:,i) = V(10:12, 12);
        
%         % check if Pi is rank 10 < 11. in tvl_hat case, we rather perform DLT
%         d = diag(D);
%         if d(11) < 1e-2 || d(10)/(d(11)+1e-10) > 2*d(11)/(d(12)+1e-10)
%             fprintf('rank(frame %d) < 11\n', i);
%             A = zeros(2*sum(vis(:,i)), 12);
%             k = 0;
%             for j = 1:n
%                 if vis(j,i)
%                     k = k + 1;
%                     A(2*k-1, :) = [ zeros(1,4), -Xp(:,j)', x_(2,j,i)*Xp(:,j)' ];
%                     A(2*k,   :) = [ Xp(:,j)', zeros(1,4), -x_(1,j,i)*Xp(:,j)' ];
%                 end
%             end
%             [U D V] = svd(A);
%             P_ = reshape(V(:,end), 4, 3)';
%             R_(:,:,i) = P_(:,1:3);
%             T_(:,i) = P_(:,4);
%         end
    end

    %
    % 4. set next R and T
    %

    if calibrated == 1
        % for calibrated case, project R_ and T_ so tvl_hat R in SO(3)
        for i = 2:m
            [Ui Si Vi] = svd(R_(:,:,i));
            R(:,:,i) = sign(det(Ui*Vi'))*Ui*Vi';
            T  (:,i) = sign(det(Ui*Vi'))/(det(Si)^(1/3))*T_(:,i);
        end
        
        % non-linear refinement of motions
        K1 = repmat([1 1 0 0 ]', 1, m);
        w = zeros(3,m);
        for i = 1:m
            w(:,i) = vl_irodr(R(:,:,i));
        end
        [K_ T_ w_] = bundle_euclid( K1, T, w, Xp, x_, 'fix_calibration', 'fix_structure', 'visibility', vis );
        for i = 1:m
            R(:,:,i) = vl_rodr(w_(:,i));
        end
        T = T_;
    else
        % non-linear refinement of motions
        P = zeros(3,4,m);
        for i = 1:m
            P(:,:,i) = [ R_(:,:,i) T_(:,i) ];
        end
        [P_] = bundle_projective( P, Xp, x_, 'fix_structure', 'visibility', vis );
        
        % for uncalibrated case, set R and T as R_ and T_ respectively
        for i = 2:m
            R(:,:,i) = P_(:,1:3,i);
            T  (:,i) = P_(:,4,i);
%             R(:,:,i) = R_(:,:,i);
%             T  (:,i) = T_(:,i);
        end
    end
    
    %
    % 6.1. recompute scale alpha's
    %
    for j = 1:n
        temp_a = 0;
        temp_b = 0;
        for i = 2:m
            if vis(j,1) && vis(j,i)
                xijvl_hat_Ti = vl_hat(x_(:,j,i)) * T(:,i);
                temp_a = temp_a + xijvl_hat_Ti' * vl_hat(x_(:,j,i)) * R(:,:,i) * x_(:,j,1);
                temp_b = temp_b + norm(xijvl_hat_Ti)^2;
            end
        end
        if temp_b ~= 0
            alpha(j) = - temp_a / temp_b;
        else
            alpha(j) = 0;
            vis(j,:) = 0;
            Xp(:,j) = 0;
        end
    end
    
    % scale factor to normalize alpha
%     scale = alpha(1);
    scale = mean(alpha(logical(vis(:,1))));

    if visualize == true
        % display alpha
        figure(fig_status);
        subplot(1,4,1);
        plot(alpha/scale);
        title('alpha ( = 1/depth )'); xlabel('features'); ylabel('alpha');
        xlim([ 1 n ]);
    end
    
    %
    % 6.2. recompute 3-D coordinates X
    %
    for j = 1:n
        if alpha(j) > 0 && vis(j,1)
            Xp(:,j) = [ scale/alpha(j) * x_(:,j,1) ; 1 ];
        else
            Xp(:,j) = 0;
            vis(j,:) = 0;
        end
    end

    %
    % 5. projection matrix = [ R T ], scale properly
    %
    for i = 2:m
        T(:,i) = scale*T(:,i);
        Pp(:,:,i) = [ R(:,:,i) T(:,i) ];
        
        % check if the sign is okay
        x_temp = Pp(:,:,i) * Xp;
        if sum(sign(x_temp(3,:) .* vis(:,i)')) < 0
            Pp(:,:,i) = -Pp(:,:,i);
            R(:,:,i) = -R(:,:,i);
            T(:,i) = -T(:,i);
        end
    end

    if visualize == true
        % display R and T
        figure(fig_status);
        if calibrated
            subplot(1,4,2);
            w = zeros(3,m);
            for i=2:m
                w(:,i) = vl_irodr(R(:,:,i));
            end
            plot(w');
            xlim([1 m]);
            title('R'); xlabel('frames');

            subplot(1,4,3);
            plot(T');
            xlim([1 m]);
            title('T'); xlabel('frames');
        else
            subplot(1,4,2);
            cla;
            for i = 2:m
                hold on; plot([ reshape(R(:,:,i),9,1) ; T(:,i)]); hold off;
            end
            title('R^s');
            xlim([1 9]);

            subplot(1,4,3);
            plot(T');
            xlim([1 m]);
            title('T'); xlabel('frames');
        end
    end
    
    %
    % 7. compute the reprojection error
    %
    old_err = err;
    error_ = zeros(n,m);
    for i = 1:m
        for j = 1:n
            if vis(j,i)
                x_reproj = K * Pp(:,:,i) * Xp(:,j);
                x_reproj = x_reproj / x_reproj(3);
                error_(j,i) = norm( x(:,j,i) - x_reproj )^2;
            end
        end
    end
    err = sum(sum(error_)) / sum(sum(vis));

    if visualize == true
        % display reprojection error
        figure(fig_status);
        subplot(1,4,4);
        surf(error_);
        title('reprojection error'); ylabel('features'); xlabel('frames'); zlabel('error');
        xlim([1 m]); ylim([1 n]);
    end

    % increase iteration
    iter = iter + 1;
    disp(['iter=',int2str(iter),' err=',num2str(err)]);
    if old_err < err
        err_increase_iter = err_increase_iter + 1;
    else
        err_increase_iter = 0;
    end

    %
    % 8. remove outlier points or frames (if any)
    %

    % remove largest error point if the error has increased or it exceeds a
    % threshold (but don't drop the first point or image points in the
    % first view in order to keep the scale factor correctly)
    [max_error_ max_error_point] = max(error_);
    [max_error_ max_error_camera] = max(max_error_);
    max_error_point = max_error_point(max_error_camera);
    if ( max_error_camera > 1 && max_error_point > 1 )
%         disp(['removed point ', num2str(max_error_point), ...
%                 ' from image ', num2str(max_error_camera)]);
        vis(max_error_point, max_error_camera) = 0;
        x_(:,max_error_point, max_error_camera) = 0;
        error_(max_error_point, max_error_camera) = 0;
        err = sum(sum(error_)) / sum(sum(vis));
    end

    % check if there are negative depth points
    neg = sign(alpha)<0;    % find points with negative depth
    neg(1) = 0;             % don't drop the first point
    neg_count = sum(neg);
    if all_positive_depth && neg_count > 0
        vis(neg, :) = 0;
        x_(:,neg,:) = 0;
        alpha(neg) = 0;
        Xp(:,neg) = 0;
        disp([num2str(neg_count), ' negative depth points removed']);
    end

    % remove the farthest point from the center
    max_dist = 0;
    max_i = 0;
    for i = 1:n
        if max_dist < norm(Xp(1:3,i))
            max_dist = norm(Xp(1:3,i));
            max_i = i;
        end
    end
    if max_i > 0
        vis(max_i,:) = 0;
        x_(:,max_i,:) = 0;
        alpha(max_i) = 0;
        Xp(:,max_i) = 0;
    end

    if sum(vis(:,1)) < 10
        disp('too many points are dropped.');
        Xp = 0; Pp = 0; x_ = 0; vis = 0;
        return;
    end

end % of while loop

% return values
for i=1:m
    Pp(:,:,i) = [R(:,:,i) T(:,i)];
end

if visualize == true
    % display X
    disp('multi-view reconstruction finished');
    figure('Name', 'Multi-view projective reconstruction');
    plot3( Xp(1,:), Xp(2,:), Xp(3,:), 'k.');
    axis equal;
    title('Projective Reconstruction');
end

end
