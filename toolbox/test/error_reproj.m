function [err error] = error_reproj( x, K, T, w, X, varargin )
% ERROR_REPROJ computes the reprojection error given K, w, T and X
%
% Input:
%       x = measurements (3-by-n-by-m)
%       K = calibration matrix or parameters (3-by-3 or 3-by-3-by-m or 4-by-m )
%       T = translation parameters (3-by-m)
%       w = exponential map parameters (3-by-m)
%       X = points in 3D (4-by-n)
%
% Options:
%        'visibility', vis : visibility map (nxm)
%
% Output:
%        err   = scalar value of the average reprojection error
%        error = n-by-m matrix of the reprojection error

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

% help
num_required_parameters = 5;
if nargin < num_required_parameters
    help error_reproj.m
    return;
end

% get information
m = size(T, 2);
n = size(X, 2);

% initialize default parameters
vis = ones(n, m);

% parse parameters
if nargin > num_required_parameters
    iVarargin = 1;
    while iVarargin <= nargin - num_required_parameters
        switch lower(varargin{iVarargin})
            case 'visibility'
                vis = varargin{iVarargin+1};
                iVarargin = iVarargin + 1;
        end
        iVarargin = iVarargin + 1;
    end
end

% get K
if size(K,1) == 4
    % make K to be 3-by-3-by-m if K is given as parameter vector
    K_ = zeros(3,3,m);
    for i = 1:m
        K_(:,:,i) = calibration_matrix( K(:,i) );
    end
    K = K_;
    clear K_;
else
    % make K to be 3-by-3-by-m if K is given as a 3-by-3 matrix
    if size(size(K), 2) == 2
        K_ = zeros(3,3,m);
        for i = 1:m
            K_(:,:,i) = K;
        end
        K = K_;
        clear K_;
    end
end

% get reprojection error
error = zeros(n, m);
for j = 1:m
    P_j = K(:,:,j) * [ vl_rodr(w(:,j)) T(:,j) ];
    for i = 1:n
        if X(4,i) == 1 && vis(n,m)
            x_reproj = P_j * X(:,i);
            x_reproj = x_reproj / x_reproj(3);
            error(i, j) = norm( x(1:2,i,j) - x_reproj(1:2) );
        end
    end
end
%error = error .* vis;
err = sum(sum(error)) / sum(sum(vis));

end
