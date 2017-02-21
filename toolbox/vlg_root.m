function path = vlg_root
% VLG_ROOT  gets VLG package root directory
%  PATH = VLGROOT() returns the root directory of the VLG
%  package.

% edited from vlfeat_root.m

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

[a,b,c] = fileparts(which('vlg_root')) ;
[a,b,c] = fileparts(a) ;
path = a;
