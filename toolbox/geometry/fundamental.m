function [ F ] = fundamental( x1, x2, varargin )
% FUNDAMENTAL calculates the Fundamental matrix from point correspondences.
%
% Input:
%        x1, x2 - image points of two views (3xn)
%
% Options:
%        'calibrated', indicates the camera is calibrated
%
% Output:
%        F - Fundamental matrix (or Essential if calibrated) (3x3)

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

% help
if nargin < 2
    help fundamental.m;
    return;
end

% get information
n = size(x1, 2);

% check minimum number of points
if n < 8
    disp('number of points should be greater or equal to 8');
    F = eye(3);
    return;
end

% initialize default parameters
calibrated = 0;

% parse parameters
if nargin > 2
    iVarargin = 1;
    while iVarargin <= nargin - 2
        switch lower(varargin{iVarargin})
            case 'calibrated'
                calibrated = 1;
        end
        iVarargin = iVarargin + 1;
    end
end

% normalize x1
x1_ = mean(x1, 2);
m_x1 = x1_(1);
m_y1 = x1_(2);
s_x1 = std(x1(1,:));
s_y1 = std(x1(2,:));
H1 = [ 1/s_x1 0 -m_x1/s_x1 ; 0 1/s_y1 -m_y1/s_y1 ; 0 0 1 ];
x1 = H1 * x1;

% normalize x2
x2_ = mean(x2, 2);
m_x2 = x2_(1);
m_y2 = x2_(2);
s_x2 = std(x2(1,:));
s_y2 = std(x2(2,:));
H2 = [ 1/s_x2 0 -m_x2/s_x2 ; 0 1/s_y2 -m_y2/s_y2 ; 0 0 1 ];
x2 = H2 * x2;

% build X
X = zeros( n, 9 );
for i = 1:n
    X(i,:) = kron(x1(:,i), x2(:,i));
end

% SVD on X = Ux Dx Vx'
[Ux Dx Vx] = svd(X);

% Fs = ninth column of Vx
Fs = Vx(:,9);

% unstack F
F_ = reshape(Fs, 3, 3);

[U D V] = svd(F_);
if calibrated == 1
    % Impose the rank=2 and s1 = s2 constraint (for calibrated)
    s = (D(1,1)+D(2,2))/2;
    D = diag([ s s 0 ]);
else
    % Impose the rank=2 constraint (for uncalibrated)
    D(3,3) = 0;
end
F_ = U * D * V';

% apply normalization of x1 and x2
F = H2' * F_ * H1;

end
