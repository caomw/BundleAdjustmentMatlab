function [Pp_ Xp_ error_] = bundle_projective( Pp, Xp, x, varargin )
% BUNDLE_PROJECTIVE performs the projective bundle adjustment.
%
% Input:
%        Po - projection matrices (3x4xm)
%        Xp - 3d points           (4xn)   in homogeneous form
%        x  - image points        (3xnxm) in homogeneous form
%
% Options:
%        'fix_structure'       - keep structure parameters fixed
%        'fix_motion'          - keep motion parameters fixed
%        'visibility', visible - set the visibility map (visible: nxm)
%        'verbose'             - display verbose information
%
% Output:
%        Pp_    - projection matrices   (3x4xm)
%        Xp_    - 3d points             (4xn)   in homogeneous form
%        error_ - errors for iterations (1x#iter)

% Copyright 2008 (C) Taehee Lee
%
% This program is part of VLG, available in the terms of the GNU
% General Public Licenseversion 2.

% help
num_required_parameters = 3;
if nargin < num_required_parameters
    help bundle_projective.m
    return;
end

% get information
m = size(Pp, 3);
n = size(x, 2);

% initialize default parameters
fix_structure   = false;
fix_motion      = false;
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
visible = double(visible); % NOTE: the mex file considers visible as double
num_vis = sum(sum(visible));

%
% build vectors: P = parameter vector in R^M, X = measurements in R^N
%

% a vector
a = zeros(12, m);
for j=1:m
    a(1:12,j) = reshape(Pp(:,:,j),12,1);
end

% b vector
b = Xp(1:3,:);

% X vector
X = x(1:2,:,:);

% assume Sx: covariance of X, to be Identity

% (i)
% init lambda
%
lambda = 0.001;

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

% while (iter < 20) && (iter<5 || error_(iter-1) - error_(iter) > 1e-3)
% while (iter < 20) && (iter<5 || error_(iter-1) - error_(iter) > 0.01*error_(iter-1))

    % initialize varaibles
    Y = zeros(12,3,n,m);

    % (ii)
    % compute derivative matrices A_ij and B_ij, and the error vector e_ij
    %
    % (iii)
    % compute intermediate expressions
    % U_j = sum_i( A_ij' * S_ij^-1 * A_ij )
    % V_i = sum_j( B_ij' * S_ij^-1 * B_ij )
    % Wij = A_ij' * S_ij^-1 * B_ij
    % eA_j = sum_i( A_ij' * S_ij^-1 * e_ij )
    % eB_i = sum_j( B_ij' * S_ij^-1 * e_ij )
    %
    [X_hat A B e U V W eA eB] = mex_bundle_proj_1_XABeUVWeAeB(a, b, X, visible);
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
        for k=1:12
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
    [S e_] = mex_bundle_proj_2_Se_(Y, W, U_, eA, eB);
    da = pinv(S) * e_;
    
    % clear memory
    clear S;

    % (vii)
    % find db by back-substitution: db_i = V*_i^-1 (eB_i - sum(W_ij'da_j) )
    %
    % (viii)
    % update the parameter vector by adding the incremental vector ( da ; db )
    % and compute the new error vector
    %
    [db a_new b_new X_hat_new] = mex_bundle_proj_3_db_new( W, da, eB, V_inv, a, b, X, visible );
    e_new = X - X_hat_new;

    e_stack = reshape(e, 2*n*m, 1);
    e_new_stack = reshape(e_new, 2*n*m, 1);
    old_error = 1/num_vis * e_stack' * e_stack;
    new_error = 1/num_vis * e_new_stack' * e_new_stack;    
    % (ix)
    % if the new error is less than the old error, then accept the new values
    % of the parameters, diminish the value of lambda by a factor of 10, and
    % start again at step (ii), or else terminate
    if new_error < old_error
        if verbose
            disp(['iter ', int2str(iter), ': error= ', ...
                  num2str(old_error), ' -> ', num2str(new_error)]);
        end
        a = a_new;
        b = b_new;
        lambda = lambda / 10;
        error_(iter) = old_error; %#ok<AGROW>
        iter = iter + 1;
        error_(iter) = new_error; %#ok<AGROW>
        iter2 = 0;
    else
    % (x)
    % if the new error is greater than the old error, then revert to the old
    % parameter values, increase the value of lambda by a factor of 10, and try
    % again from step (iv)
        lambda = lambda * 10;
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

% projection
Pp_ = zeros(3,4,m);
for j=1:m
    Pp_(:,:,j) = reshape(a(:,j),3,4);
end
% structure
Xp_ = [ b ; Xp(4,:) ];

end
