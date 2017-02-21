function [K_ Te_ w_ Xe_ error_] = bundle_euclid_nomex( K, Te, w, Xe, x, varargin )
% BUNDLE_EUCLID_NOMEX  performs Euclidean bundle adjustment (without MEX)
%
% Input:
%        K  - calibration parameters (4xm) ([fx fy cx cy]' x m-images)
%        Te - translation            (3xm)
%        w  - rotation               (3xm)
%        Xe - 3d points              (4xn)   in homogeneous form
%        x  - image points           (3xnxm) in homoegeneous form
%
% Options:
%        'fix_structure'       - keep structure parameters fixed
%        'fix_motion'          - keep motion parameters fixed
%        'fix_calibration'     - keep calibration parameters fixed
%        'fix_principal'       - keep principal points fixed
%        'visibility', visible - set the visibility map (visible: nxm)
%        'verbose'             - display verbose information
%
% Output:
%        K_     - calibration parameters (4xm)
%        Te_    - translation            (3xm)
%        w_     - rotation               (3xm)
%        Xe_    - 3d points              (4xn) in homogeneous form
%        error_ - errors for iterations  (1x#iter)

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

% help
num_required_parameters = 5;
if nargin < num_required_parameters
    help bundle_euclid_nomex.m
    return;
end

% get information
m = size(w, 2);
n = size(x, 2);

% initialize default parameters
fix_structure   = false;
fix_motion      = false;
num_variableK   = 4;
visible         = reshape(x(1,:,:)~=0 | x(2,:,:)~=0, n, m);
verbose         = false;

% parse optional parameters
if nargin > num_required_parameters
    iVarargin = 1;
    while iVarargin <= nargin - num_required_parameters
        switch lower(varargin{iVarargin})
            case 'fix_structure'
                fix_structure = true;
            case 'fix_motion'
                fix_motion = true;
            case 'fix_calibration'
                num_variableK = 0;
            case 'fix_principal'
                num_variableK = 1;
            case 'visibility'
                visible = varargin{iVarargin+1};
                iVarargin = iVarargin + 1;
            case 'verbose'
                verbose = true;
        end
        iVarargin = iVarargin + 1;
    end
end

% count the visible measurements
num_vis = sum(sum(visible));

% a vector
num_a = 6+num_variableK;
a = zeros(num_a, m);
a(1:6,:) = [ w ; Te ];
if num_variableK == 1
    a(7,:) = K(1,:);
elseif num_variableK == 4
    a(7:end,:) = K(:,:);
end

% b vector
b = Xe(1:3,:);

% X vector
X = x(1:2,:,:);

% assume Sx: covariance of X, to be Identity


% (i)
% init lambda
%
lambda = 0.001;
nu = 2;

% while loop
iter = 1;
iter2 = 0;
max_iter = 20;
max_iter2 = 10;
error_ = [];
while (iter  < max_iter ) && ...
      (iter2 < max_iter2) && ...
      (iter  < 3 || ...
      (error_(iter)>1e-20 && error_(iter-1)-error_(iter) > 1e-3*error_(iter-1)))

    % init variables
    X_hat = zeros(2, n, m);
    X_hat_new = zeros(2, n, m);
    e = zeros(2, n, m);
    A = zeros(2,num_a,n,m);
    B = zeros(2,3,n,m);
    Y = zeros(num_a,3,n,m);
   
    % (ii)
    % compute derivative matrices
    %     A_ij and B_ij, and the error vector e_ij
    %
    for i=1:n           % for each point i
        for j=1:m       % for each camera j
            if visible(i,j)
                % reprojection of point
                X_hat(:,i,j) = reprojection_point( ...
                    K(:,j), a(:,j), b(:,i), num_variableK );

                % A_ij
                for k=1:num_a
                    dir = zeros(num_a, 1);
                    dir(k,1) = 1;
                    A(:,k,i,j) = derivative_camera( ...
                        a(:,j), dir, X_hat(:,i,j), K(:,j), b(:,i), num_variableK );
                end
                
                % B_ij
                for k=1:3
                    dir = zeros(3, 1);
                    dir(k,1) = 1;
                    B(:,k,i,j) = derivative_point( ...
                        b(:,i), dir, X_hat(:,i,j), K(:,j), a(:,j), num_variableK );
                end

                % e_ij
                e(:,i,j) = X(:,i,j) - X_hat(:,i,j);
            else
                X_hat(:,i,j) = X(:,i,j);
                A(:,:,i,j) = 0;
                B(:,:,i,j) = 0;
                e(:,i,j)   = 0;
            end
        end
    end
    
    % (iii)
    % compute intermediate expressions
    %     U_j = sum_i( A_ij' * S_ij^-1 * A_ij )
    %     V_i = sum_j( B_ij' * S_ij^-1 * B_ij )
    %     Wij = A_ij' * S_ij^-1 * B_ij
    %     eA_j = sum_i( A_ij' * S_ij^-1 * e_ij )
    %     eB_i = sum_j( B_ij' * S_ij^-1 * e_ij )
    %
    U = zeros(num_a,num_a,m);
    V = zeros(3,3,n);
    W = zeros(num_a,3,n,m);
    eA = zeros(num_a,m);
    eB = zeros(3,n);
    for i=1:n
        for j=1:m
            U(:,:,j) = U(:,:,j) + A(:,:,i,j)' * A(:,:,i,j);
            V(:,:,i) = V(:,:,i) + B(:,:,i,j)' * B(:,:,i,j);
            W(:,:,i,j) = A(:,:,i,j)' * B(:,:,i,j);
            eA(:,j) = eA(:,j) + A(:,:,i,j)' * e(:,i,j);
            eB(:,i) = eB(:,i) + B(:,:,i,j)' * e(:,i,j);
        end
    end

    if fix_structure
        V (:,:,:)   = 0;
        W (:,:,:,:) = 0;
        eB(:,:)     = 0;
    end
    if fix_motion
        U (:,:,:)   = 0;
        W (:,:,:,:) = 0;
        eA(:,:)     = 0;
    end
    
    % clear memory
    clear A B;

    % (iv)
    % augment U and V by multiplying their diagonal elements by 1+lambda
    %
    U_ = U;
    for j=1:m
        for k=1:num_a
            U_(k,k,j) = (1 + lambda) * U(k,k,j);
        end
    end
    V_ = V;
    for i=1:n
        for k=1:3
            V_(k,k,i) = (1 + lambda) * V(k,k,i);
        end
    end

    % (v)
    % Y_ij = W_ij V_i*^-1
    %
    V_inv = zeros(3,3,n);
    for i=1:n
        V_inv(:,:,i) = pinv(V_(:,:,i));
        for j=1:m
            Y(:,:,i,j) = W(:,:,i,j) * V_inv(:,:,i);
        end
    end

    % clear memory
    clear U V V_;
    
    % (vi)
    % find da by solving S da = e_
    %
    S = zeros(num_a*m, num_a*m);
    for j=1:m
        for k=1:m
            Sjk = zeros(num_a,num_a);
            if j == k
                % S_jj = - sum_i( Y_ij W_ij' + U*_j ) + U*_(j)
                Sjk = U_(:,:,j);
                for i=1:n
                    Sjk = Sjk - Y(:,:,i,j) * W(:,:,i,j)';
                end
            else
                % S_jk = - sum_i( Y_ij W_ik' )
                for i=1:n
                    Sjk = Sjk - Y(:,:,i,j) * W(:,:,i,k)';
                end
            end
            j1=num_a*(j-1)+1;
            k1=num_a*(k-1)+1;
            S(j1:(j1+num_a-1), k1:(k1+num_a-1)) = Sjk;
        end
    end
    e_ = zeros(num_a*m, 1);
    k = 1;
    for j=1:m
        Sum_Y_eB = zeros(num_a,1);
        for i=1:n
            Sum_Y_eB = Sum_Y_eB + Y(:,:,i,j) * eB(:,i);
        end
        e_(k:(k+num_a-1),1) = eA(:,j) - Sum_Y_eB;
        k = k + num_a;
    end
    da = pinv(S) * e_;

    % clear memory
    clear S;
    
    % (vii)
    % find db by back-substitution: db_i = V*_i^-1 (eB_i - sum(W_ij'da_j) )
    % (viii)
    % update the parameter vector by adding the incremental vector ( da ; db )
    % and compute the new error vector
    %
    db = zeros(3,n);
    for i=1:n
        Sum_W_da = zeros(3,1);
        k = 1;
        for j=1:m
            Sum_W_da = Sum_W_da + W(:,:,i,j)' * da(k:(k+num_a-1), 1);
            k = k + num_a;
        end
        db(:,i) = V_inv(:,:,i) * ( eB(:,i) - Sum_W_da );
    end

    % (viii)
    % update the parameter vector by adding the incremental vector ( da ; db )
    % and compute the new error vector
    %
    a_new = zeros(num_a,m);
    k = 1;
    for j=1:m
        a_new(:,j) = a(:,j) + da(k:(k+num_a-1),1);
        k = k+num_a;
    end
    b_new = b + db;
    for i=1:n
        for j=1:m
            % new reprojection of point
            if visible(i,j)
                X_hat_new(:,i,j) = reprojection_point( ...
                    K(:,j), a_new(:,j), b_new(:,i), num_variableK );
            else
                X_hat_new(:,i,j) = X(:,i,j);
            end
        end
    end
    e_new = X - X_hat_new;

    e_stack = reshape(e, 2*n*m, 1);
    e_new_stack = reshape(e_new, 2*n*m, 1);
    old_error = e_stack' * e_stack;
    new_error = e_new_stack' * e_new_stack;    
    % (ix)
    % if the new error is less than the old error, then accept the new values
    % of the parameters, diminish the value of lambda by a factor of 10, and
    % start again at step (ii), or else terminate
    g = [ reshape( eA , num_a*m, 1 ) ; reshape( eB , 3*n, 1 ) ];
    dp = [ da ; reshape(db, 3*n, 1) ];
    rho = ( old_error - new_error ) / ( dp' * ( lambda * dp + g ) );
    if ( old_error - new_error ) > 0
        old_error = old_error / num_vis;
        new_error = new_error / num_vis;
        if verbose
            disp(['iter ', int2str(iter), ': error= ', ...
                  num2str(old_error), ' -> ', num2str(new_error)]);
        end
        a = a_new;
        b = b_new;
        lambda = lambda * max( 1/3, 1-(2*rho-1)^3);
        nu = 2;
        error_(iter) = old_error; %#ok<AGROW>
        iter = iter + 1;
        error_(iter) = new_error; %#ok<AGROW>
        iter2 = 0;
    else
    % (x)
    % if the new error is greater than the old error, then revert to the old
    % parameter values, increase the value of lambda by a factor of 10, and try
    % again from step (iv)
        lambda = lambda * nu;
        nu = 2 * nu;
        iter2 = iter2 + 1;
    end
    
    % clear memory
    clear X_hat A B e U V W eA eB;
    clear U_ V_ V_inv Y;
    clear S e_ da db a_new b_new X_hat_new e_new;
    clear e_stack e_new_stack old_error new_error;

end % of while loop

%
% build output result
%

% calibration matrix
K_  = K;
if num_variableK == 1
    K_(1,:) = a(7,:);
    K_(2,:) = a(7,:);
elseif num_variableK == 4
    K_(:,:) = a(7:end,:);
end
% motion
w_  = a(1:3,:);
Te_ = a(4:6,:);
% structure
Xe_ = [ b ; ones(1,n) ];

end
