function H = display_image_pair( I1, I2, x1, x2 )
% DISPLAY_IMAGE_PAIR displays a pair of image with their feature correspondences.
%
% I1, I2: images
% x1, x2: 3-by-n for n features
%
% return the handle of the figure, H

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

index = x1(1,:) ~= 0 & x1(2,:) ~= 0 & x2(1,:) ~= 0 & x2(2,:) ~= 0;
x1 = x1(:,index);
x2 = x2(:,index);

% get information
n = size(x1, 2);
w = size(I1, 2);

I = [ I1 I2 ];

H = figure('Name', 'Image pair with point correspondences');

imshow(I);
hold on;
axis off; axis ij; axis equal;
for i = 1:n
    line([x1(1,i) x2(1,i)+size(I1,2)], [x1(2,i) x2(2,i)]);
end
plot(x1(1,:), x1(2,:), 'g.');
plot(x2(1,:)+w, x2(2,:), 'g.');

hold off;
title('Image pair with point correspondences');

end
