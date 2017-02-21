function H = display_image_pair_with_F( I1, I2, F, x1, x2, varargin )
% DISPLAY_IMAGE_PAIR_WITH_F displays a pair of image with their features and the Fundamental matrix.
%
% I1, I2 : images
% F      : Fundamental matrix I1 -> I2
% x1, x2 : features (3-by-n for n features)
%
% (optional)
% range  : range of image plane ( 2x2 for [ left right ; bottom top ] )
%
% return the figure's handle H
%

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

% get information
n = size(x1, 2);
range = ones(2,2);
if nargin == 6
    range = varargin{1};
else
    range = [ 1 size(I1, 2) ; 1 size(I1, 1) ];
end
w = range(1,2) - range(1,1);

I = [ I1 I2 ];

H = figure('Name', 'Image pair with Epipolar lines');

imshow(I); axis ij; axis equal;
hold on;
for i = 1:n
    display_line( F'*x2(:,i), range(1,:), range(2,:)          );
    display_line( F*x1(:,i),  range(1,:), range(2,:), [ w 0 ] );
end
plot(x1(1,:),   x1(2,:), 'g.');
plot(x2(1,:)+w, x2(2,:), 'g.');

hold off;
title('Image pair with Epipolar lines');

% subplot(1,2,1); axis ij; axis equal;
% imshow(I1);
% hold on;
% for i = 1:n
%     display_line( F'*x2(:,i),  range(1,:), range(2,:) );
% end
% plot(x1(1,:), x1(2,:), 'g.');
% hold off;
% 
% subplot(1,2,2); axis ij; axis equal;
% imshow(I2);
% hold on;
% for i = 1:n
%     display_line( F*x1(:,i), range(1,:), range(2,:) );
% end
% plot(x2(1,:), x2(2,:), 'g.');
% hold off;

end
