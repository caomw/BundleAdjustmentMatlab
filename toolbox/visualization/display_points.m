function display_points( X, varargin )
% DISPLAY_POINTS displays 3D points.
%
% Input:
%       X - 3d location of points (3xn) or (4xn)
%
% Options:
%       'pixelsize', pixelsize - pixel size
%       'color', color         - colors nx3 in RGB vector

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

% help
if nargin < 1
    help display_points.m
    return;
end

% get information
n = size(X,2);

% initialize default parameters
pixelsize = ones(1,n);
color = zeros(n,3);

% parse parameters
if nargin > 1
    iVarargin = 1;
    while iVarargin <= nargin - 1
        switch lower(varargin{iVarargin})
            case 'pixelsize'
                pixelsize = varargin{iVarargin+1}*ones(1,n);
                iVarargin = iVarargin + 1;
            case 'color'
                color = varargin{iVarargin+1};
                if size(color,1) == 1
                    color = repmat(color, n,1);
                end
                iVarargin = iVarargin + 1;
        end
        iVarargin = iVarargin + 1;
    end
end

% select only valid points if X is given as 4xn form
if size(X,1) == 4
    index = X(4,:)==1;
    X = X(1:3,index);
    pixelsize = pixelsize(:,index);
    color = color(index,:);
end

% plot points
hold on;
scatter3(X(1,:), X(2,:), X(3,:), pixelsize, color, 'filled'), axis equal;
hold off;

end
