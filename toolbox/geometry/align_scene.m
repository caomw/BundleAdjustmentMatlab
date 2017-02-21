function [ T_ Omega_ X_ ] = align_scene( T, Omega, X, TRef, OmegaRef, XRef, ScaleOption )
% ALIGN_SCENE aligns a scene.
%
% This aligns a scene represented by (T, Omega, X) to a reference frame
% (TRef, OmegaRef, XRef). The rotation and translation is aligned so that
% the first frame of T and Omega is same as the first frame of TRef and
% OmegaRef. The scale factor is aligned so that the norm of the norm of
% centroid of X is same as the norm of centroid of XRef.
%
% When TRef, OmegaRef and XRef are omitted, by default, reference frame
% becomes [R|T] = [I|0], and the norm of centroid of 3d points = 1.
% For example,
% >> [T_ Omega_ X_] = align_scene( T, Omega, X );
%
% Input:
%        T        - camera translation           (3xm)
%        Omega    - camera rotation              (3xm)
%        X        - 3D point location            (4xn) in homogeneous form
% (Optional)
%        TRef     - reference camera translation (3xm)
%        OmegaRef - reference camera rotation    (3xm)
%        XRef     - reference 3d point location  (4xn) in homogeneous form
%
%        ScaleOption - set the scale factor ( 'centroid', 'translation' )
%
% Output:
%        T_       - aligned camera translation   (3xm)
%        Omega_   - aligned camera rotation      (3xm)
%        X_       - aligned 3D point location    (4xn) in homogeneous form

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

% help
if nargin ~= 3 && nargin ~= 6 && nargin ~= 7
    help align_scene.m;
    return;
end

% get information
m = size(T,2);
n = size(X,2);
if nargin ~= 7
    ScaleOption = 'centroid';
    switch lower(ScaleOption)
        case 'centroid'
            ScaleOption = 'centroid';
        case 'translation'
            ScaleOption = 'translation';
        otherwise
            ScaleOption = 'centroid';
    end
end
iRef = 2;

% compute the reference frame's transformation
if nargin == 3
    TRef     = zeros(3,m);
    OmegaRef = zeros(3,m);
    XRef     = ones(4,n);
end
Tref1 = TRef(:,1);
Rref1 = vl_rodr( OmegaRef(:,1) );
Xref1 = [ (Rref1*XRef(1:3,:)+repmat(Tref1,1,size(XRef,2))).*repmat(XRef(4,:),3,1) ; XRef(4,:) ];
switch ScaleOption
    case 'centroid'
        Sref  = 1/norm(mean(Xref1(1:3,XRef(4,:)==1),2));
    case 'translation'
        Sref  = 1/norm(Rref1'*Tref1 - vl_rodr(OmegaRef(:,iRef))'*TRef(:,iRef));
end

% align the rotation and translation
T_     = zeros(3, m);
Omega_ = zeros(3, m);
R1 = vl_rodr(Omega(:,1));
if isnan(R1) | isinf(R1)
    keyboard;
end

T1 = T(:,1);
X1 = [ (R1*X(1:3,:)+repmat(T1,1,n)).*repmat(X(4,:),3,1) ; X(4,:) ];
switch ScaleOption
    case 'centroid'
        S1 = 1/norm(mean(X1(1:3,X(4,:)==1),2));
    case 'translation'
        S1 = 1/norm(R1'*T1 - vl_rodr(Omega(:,iRef))'*T(:,iRef));
end
for j = 1:m
    RjR1_ = vl_rodr(Omega(:,j)) * R1';
    Omega_(:,j) = vl_irodr( RjR1_ * Rref1 );
    if isnan(Omega_(:,j))
        Omega_(:,j) = rodrigues( RjR1_ * Rref1 );
    end
    T_(:,j) = RjR1_*Tref1 + S1/Sref*( - RjR1_ * T1 + T(:,j) );
end
X_ = zeros(4, n);
for i = 1:n
    if X(4,i) == 1
        X_(1:3,i) = Rref1'*(S1/Sref*(R1*X(1:3,i)+T1) - Tref1);
        X_(  4,i) = 1;
    end
end

end
