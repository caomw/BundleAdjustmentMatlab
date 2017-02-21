function [ Xe Re Te K ] = upgrade_euclid( Xp, Pp )
% UPGRADE_EUCLID upgrades to the euclidean reconstruction
%
% This function upgrades the projective reconstruction to the euclidean
% reconstruction by recovering the absolute quadric.
%
% Input:
%       Xp - projective 3d location of points (4xn)
%       Pp - projection matrices              (3x4xm)
%
% Output:
%       Xe - Euclidean 3d location of points (4xn)
%       Re - camera rotation matrix          (3x3xm)
%       Te - camera translation vector       (3xm)
%       K  - camera calibration matrix       (3x3xm)

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

% help
if nargin ~= 2
    help upgrade_euclid.m;
    return;
end

% get information
m = size(Pp, 3);
n = size(Xp, 2);

%
% 2. form Chi (4m-by-5)
% 3. form b   (4m)
%
Chi = [];
b = [];
for i = 1:m
    u = Pp(1,:,i); u1 = u(1); u2 = u(2); u3 = u(3); u4 = u(4);
    v = Pp(2,:,i); v1 = v(1); v2 = v(2); v3 = v(3); v4 = v(4);
    w = Pp(3,:,i); w1 = w(1); w2 = w(2); w3 = w(3); w4 = w(4);

    % Chi
    Chi = [ Chi ; ...
        u1^2+u2^2-v1^2-v2^2 2*u4*u1-2*v1*v4 2*u4*u2-2*v2*v4 2*u4*u3-2*v3*v4 u4^2-v4^2 ; ...
        u1*v1+u2*v2         u4*v1+u1*v4     u4*v2+u2*v4     u4*v3+u3*v4     u4*v4     ; ...
        u1*w1+u2*w2         u4*w1+u1*w4     u4*w2+u2*w4     u4*w3+u3*w4     u4*w4     ; ...
        v1*w1+v2*w2         v4*w1+v1*w4     v4*w2+v2*w4     v4*w3+v3*w4     v4*w4       ...
        ]; %#ok<AGROW>
    
    % b
    b = [ b ; -u3^2+v3^2 ; -u3*v3 ; -u3*w3 ; -v3*w3 ]; %#ok<AGROW>
end

%
% 4. solve Qs = Chi \ b
%
Qs = Chi \ b;

%
% 5. unstack Q_
%
Q_ = [
    Qs(1) 0     0     Qs(2) ; ...
    0     Qs(1) 0     Qs(3) ; ...
    0     0     1     Qs(4) ; ...
    Qs(2) Qs(3) Qs(4) Qs(5)   ...
    ];

%
% 6. enforce the rank-3 constraint on Q_
%
[U D V] = svd(Q_);
disp(diag(D));
D(4,4) = 0;
Q = U * D * V';
if Q(1,1) < 0
    disp('warning: Q(1,1) is negative.');
    Q = - Q;
end

% %
% % 7. compute focal lengths
% %
% f2 = zeros(3,3,m);
% for i=1:m
%     f2(:,:,i) = Pp(:,:,i) * Q * Pp(:,:,i)';
%     f2(:,:,i) = f2(:,:,i) / f2(3,3,i);
% end

%
% 8.1. construct K1 and v
%
a1 = Q(1,1);
a2 = Q(1,4);
a3 = Q(2,4);
a4 = Q(3,4);
K1 = eye(3);
K1(1,1) = sqrt(a1);
K1(2,2) = sqrt(a1);
v = - [ a2/a1 ; a3/a1 ; a4 ];

%
% 8.2. euclidean upgrade with H
%
H = [ K1 zeros(3,1) ; -v'*K1 1 ];
Xe = inv(H) * Xp;
for i = 1:n
    if Xp(4,i) > 0
        Xe(:,i) = Xe(:,i) / Xe(4,i);
    end
end
% make positive depths
if ( sum(sign(Xe(3,:))<0) >= n/2 )
    disp('flipping the sign to ensure positive depths');
    % flip the H sign and compute Xe again
    H(4,4) = -H(4,4);
    Xe = inv(H) * Xp;
    for i = 1:n
        if Xp(4,i) > 0
            Xe(:,i) = Xe(:,i) / Xe(4,i);
        end
    end
    % check if any negative depth points remain.
    disp([num2str(sum(sign(Xe(3,:))>0)), ' positive depth points.']);
    if ( sum(sign(Xe(3,:))<0 ) > 0 )
        disp([num2str(sum(sign(Xe(3,:))<0)), ' negative depth point remain.']);
    end
end

% reconstruct rigid camera motion
K  = zeros(3,3,m);
Re = zeros(3,3,m);
Te = zeros(3,m);
% err = 0;

for i=1:m
    Pe = Pp(:,:,i) * H;
    
    K (:,:,i) = K1;
    Re(:,:,i) = vl_rodr(vl_irodr( K1 \ Pe(1:3,1:3) ));
    Te  (:,i) = K1 \ Pe(1:3,4);

    if sum(sum(isnan(Re(:,:,i)))) > 0 || sum(isnan(Te(:,i))) > 0
        % factorize K, R and T
        [r q] = rq(Pe(1:3,1:3));
        K(:,:,i) = r;
        Re(:,:,i) = q;
        Te(:,i) = r \ Pe(1:3,4);
        % normalize K to a calibration matrix
        K33_factor = K(3,3,i);
        K(:,:,i) = K(:,:,i) / K33_factor;

        % make the focal length positive
        if K(1,1,i) < 0
            K (:,:,i) = K(:,:,i) * diag([-1,1,1]);
            Re(:,:,i) = diag([-1,1,1]) * Re(:,:,i);
            Te  (:,i) = diag([-1,1,1]) * Te  (:,i);
        end
        if K(2,2,i) < 0
            K (:,:,i) = K(:,:,i) * diag([1,-1,1]);
            Re(:,:,i) = diag([1,-1,1]) * Re(:,:,i);
            Te  (:,i) = diag([1,-1,1]) * Te  (:,i);
        end
        % make the rotation matrix's determinant = 1
        if det(Re(:,:,i)) < 0
            Re(:,:,i) = -Re(:,:,i);
            Te  (:,i) = -Te  (:,i);
        end
        % still, use the calibration matrix from above
        K(:,:,i) = K1;
    end
    
%     for j=1:n
%         x_reproj = Kguess * K(:,:,i) * [Re(:,:,i) Te(:,i)] * Xe(:,j);
%         x_reproj = x_reproj / x_reproj(3);
%         err = err + norm(x(:,j,i) - x_reproj)^2;
%     end
end
% err = err / m / n;
% disp(['error after upgrade = ', num2str(err)]);

% % display euclidean reconstruction
% disp('euclidean upgrade finished');
% figure('Name', 'Initial Euclidean reconstruction');
% display_reconstruction( Xe, 'k.', Re, Te, 'g.', 1 );
% title('Euclidean Reconstruction (Before refinement)');

end
