function [ matches ] = match_sift_unique( d1, d2 )
% MATCH_SIFT_UNIQUE matches two sets of SIFT descriptors.
%
% Input:
%        d1, d2 - SIFT descriptors from two images
%
% Output:
%        matches - matched feature indices (2xn)

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

% ubcmatch between d1 and d2
matches_ubc = vl_ubcmatch( d1, d2 );

% remove the points that are matched to a same 3d point
n_match_ubc = size(matches_ubc,2);
matches = zeros(2,n_match_ubc);
n_match = 0;
for i = 1:n_match_ubc
    search = matches_ubc(2,:) == matches_ubc(2,i);
    if sum(search) == 1
        % if there is no duplicated points, keep it
        n_match = n_match + 1;
        matches(:,n_match) = matches_ubc(:,i);
    end
end
% cut off duplicated (thus dropped) matches
matches = matches(:,1:n_match);

end
