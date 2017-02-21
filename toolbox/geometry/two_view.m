function [ X lambda R T F error ] = two_view( x1, x2, varargin )
% TWO_VIEW performs a two-view projective reconstruction algorithm.
%
% Input:
%        x1, x2 - image points of two views (3xn)
%
% Options:
%        'calibrated' - indicates the camera is calibrated
%        'K', K       - calibration matrix (3x3)
%
% Output:
%        X      - 3d (projective) location of poitns (4xn)
%        lambda - depths of the points               (1xn)
%        R      - camera rotation                    (3x3)
%        T      - camera translation                 (3x1)
%        F      - fundamental matrix (or essential)  (3x3)
%        error  - reprojection error

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

% help
num_required_parameters = 2;
if nargin < num_required_parameters
    help two_view.m
    return;
end

% initialize default parameters
calibrated = false;
K = eye(3);

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
        end
        iVarargin = iVarargin + 1;
    end
end

% get information
n = size(x1, 2);
if n < 8
    disp('not enough point correspondences');
    X = 0; lambda = 0; R = 0 ; T = 0; F = 0; error = 0;
    return;
end

% normalize x1 and x2
Kinv = inv(K);
x1_ = Kinv * x1;
x2_ = Kinv * x2;

% compute the Fundamental matrix
if calibrated
    F_ = fundamental( x1_, x2_, 'calibrated' );
else
    F_ = fundamental( x1_, x2_ );
end

% denormalization
F = Kinv' * F_ * Kinv;

% decompose F into R and T (taking F as if it is the Essential matrix)
[U D V] = svd(F_);

% make possible solutions
Rz = [ 0 1 0 ; -1 0 0 ; 0 0 1 ];
P2(:,:,1) = [ U*Rz*V'  vl_ihat(U*Rz'*diag([1,1,0])*U') ];
P2(:,:,2) = [ U*Rz*V'  vl_ihat(U*Rz *diag([1,1,0])*U') ];
P2(:,:,3) = [ U*Rz'*V' vl_ihat(U*Rz'*diag([1,1,0])*U') ];
P2(:,:,4) = [ U*Rz'*V' vl_ihat(U*Rz *diag([1,1,0])*U') ];

[U D V] = svd(-F_);
P2(:,:,5) = [ U*Rz*V'  vl_ihat(U*Rz'*diag([1,1,0])*U') ];
P2(:,:,6) = [ U*Rz*V'  vl_ihat(U*Rz *diag([1,1,0])*U') ];
P2(:,:,7) = [ U*Rz'*V' vl_ihat(U*Rz'*diag([1,1,0])*U') ];
P2(:,:,8) = [ U*Rz'*V' vl_ihat(U*Rz *diag([1,1,0])*U') ];

% compute 3d projective structure X for each solution
X = zeros(4,n);
R = eye(3);
T = zeros(3,1);
max_positive_depths = 0;
max_k = 0;
for k = 1:8
    R_ = P2(:,1:3,k);
    T_ = P2(:,4,k);

    if det(R_) > 0.999
        for j = 1:n
            M = [ -1  0 x1_(1,j) 0 ; ...
                   0 -1 x1_(2,j) 0 ; ...
                  (x2_(1,j)*R_(3,:) - R_(1,:)) (x2_(1,j)*T_(3) - T_(1)) ; ...
                  (x2_(2,j)*R_(3,:) - R_(2,:)) (x2_(2,j)*T_(3) - T_(2)) ];
            [Up Dp Vp] = svd( M );
            Xp_(:,j) = Vp(:,4) / Vp(4,4); %#ok<AGROW>
        end
        % compute depths
        xp2 = [R_ T_] * Xp_;
        positive_depths = sum(sign(Xp_(3,:))>0) + ...
                          sum(sign(xp2(3,:))>0) ;
        if (max_positive_depths < positive_depths)
            max_k = k;
            max_positive_depths = positive_depths;
            X = Xp_;
            R = R_;
            T = T_;
        end
    end
end
lambda = X(3, :);

if max_k > 0
    % compute reprojection error
    error_ = zeros(2, n);
    for j = 1:n
        % first frame
        x_reproj = K * [ eye(3,3) zeros(3,1) ] * X(:,j);
        x_reproj = x_reproj / x_reproj(3);
        error_(1,j) = norm( x1(:,j) - x_reproj )^2;

        % second frame
        x_reproj = K * [ R T ] * X(:,j);
        x_reproj = x_reproj / x_reproj(3);
        error_(2,j) = norm( x2(:,j) - x_reproj )^2;
    end
    error = sum(sum(error_)) / 2 / n;

    disp(['(error = ',...
          num2str(sum(error_(1,:))/n), ' , ', ...
          num2str(sum(error_(2,:))/n), ' )']);
else
    error = 0;
    disp([num2str(floor(max_positive_depths/2)), ' positive depths among ',...
          num2str(n), ' in two views']);
    disp('failed to get a physically possible cameras');
end

end
