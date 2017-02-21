function display_line( l, w, h, varargin )
% DISPLAY_LINE draws a line within a rectangle defined by w and h
%
% l: 3-vector representing a line
% w: 2-vector representing x-direction range: e.g. [0 640]
% h: 2-vector representing y-direction range: e.g. [0 480]
%
% (optional)
% shift: shifting the line to [ w h ]
%

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

shift = [ 0 0 ];
if nargin  == 4
    shift = varargin{1};
end

visible = 1;
if ( l(2) ~= 0 ) && ( l(1) ~= 0 )
    % lines neither horizontal or vertical
    x = w;
    y = -(l(1)/l(2)).*x - (l(3)/l(2));
    
    % take the range into account
    if y(1) < h(1) || y(1) > h(2)
        y(1) = min(max(y(1), h(1)), h(2));
        x(1) = -(l(2)/l(1))*y(1) - (l(3)/l(1));
    end
    if y(2) < h(1) || y(2) > h(2)
        y(2) = min(max(y(2), h(1)), h(2));
        x(2) = -(l(2)/l(1))*y(2) - (l(3)/l(1));
    end
else
    if l(2) ~= 0
        % l(1) == 0, horizontal line
        x = w;
        y = [ -l(3)/l(2) -l(3)/l(2) ];
    elseif l(1) ~= 0
        % l(2) == 0, vertical line
        y = h;
        x = [ -l(3)/l(1) -l(3)/l(1) ];
    else
        visible = 0;
    end
end

% plot only visible line
if visible == 1
    plot(x+shift(1), y+shift(2));
end

end
