function path = vlg_setup
% VLG_SETUP adds VLG toolbox path to MATLAB path
%  PATH = VLG_SETUP() adds VLG to MATLAB path.

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

root=vlg_root ;
addpath(fullfile(root,'toolbox'                )) ;
addpath(fullfile(root,'toolbox','bundle'       )) ;
addpath(fullfile(root,'toolbox','geometry'     )) ;
addpath(fullfile(root,'toolbox','visualization')) ;
addpath(fullfile(root,'toolbox','test'         )) ;

fprintf('Welcome to VLG!\n') ;

