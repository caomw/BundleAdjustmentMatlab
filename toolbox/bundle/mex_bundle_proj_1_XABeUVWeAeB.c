/*
 * mex_bundle_proj_1_XABeUVWeAeB.c
 *
 * Copyright 2008 (C) Taehee Lee
 *
 * This program is part of VLG, available in the terms of the GNU
 * General Public Licenseversion 2.
 */
#include <mexutils.h>

#define MAX_NUM_A   12

void reproject_projective_point(
    const double a[],
    const double b[],
    double x[2])
{
    double P[12];
    double x_[3];

    /* Construct P */
    memcpy(P, a, sizeof(double)*12);
    
    /* x_ = P * [ b ; 1 ] */
    x_[0] = P[0]*b[0] + P[3]*b[1] + P[6]*b[2] + P[9];
    x_[1] = P[1]*b[0] + P[4]*b[1] + P[7]*b[2] + P[10];
    x_[2] = P[2]*b[0] + P[5]*b[1] + P[8]*b[2] + P[11];

    /* x = x_(1:2) / x_(3) */
    x[0] = x_[0] / x_[2];
    x[1] = x_[1] / x_[2];
}

void derivative_projective_camera(
    const double a0[12],
    const double da[12],
    const double X0[2],
    const double b[3],
    double dXda[2])
{
    double h = 1e-10;
    double a1[12];
    double X1[2];
    
    int i;

    /* a1 = a0 + h * da */
    for ( i = 0 ; i < 12 ; i ++ )
    {
        a1[i] = a0[i] + h * da[i];
    }
    
    /* f1 = reprojection_projective_point( a1, b ) */
    reproject_projective_point( a1, b, X1 );
    
    /* dXda = (f1 - f0) / h */
    dXda[0] = (X1[0] - X0[0]) / h;
    dXda[1] = (X1[1] - X0[1]) / h;
}

void derivative_projective_point(
    const double b0[3],
    const double db[3],
    const double X0[2],
    const double a[12],
    double dXdb[2])
{
    double h = 1e-10;
    double b1[3];
    double X1[2];
    
    int i;

    /* b1 = b0 + h * db */
    for ( i = 0 ; i < 3 ; i ++ )
    {
        b1[i] = b0[i] + h * db[i];
    }
    
    /* f1 = reprojection_projective_point( a, b1 ) */
    reproject_projective_point( a, b1, X1 );
    
    /* dXdb = (f1 - f0) / h */
    dXdb[0] = (X1[0] - X0[0]) / h;
    dXdb[1] = (X1[1] - X0[1]) / h;
}

void mexFunction( int nout, mxArray *pout[],
                  int nin,  const mxArray *pin[] )
{
    int i, j, k, row, col;

    /*
     * get parameters
     * a: pin[0]        (12xm)
     * b: pin[1]        (3xn)
     * X: pin[2]        (2xnxm)
     * visible: pin[3]  (nxm)
     */
    double *pa = mxGetPr(pin[0]);
    double *pb = mxGetPr(pin[1]);
    double *pX = mxGetPr(pin[2]);
    double *pvisible = mxGetPr(pin[3]);
    int m, n;
    int num_a;

    /*
     * output parameters
     */
    mwSize dim_X_hat[3];
    mwSize dim_A[4];
    mwSize dim_B[4];
    mwSize dim_e[3];
    mwSize dim_U[3];
    mwSize dim_V[3];
    mwSize dim_W[4];
    double *X_hat;
    double *A;
    double *B;
    double *e;
    double *U;
    double *V;
    double *W;
    double *eA;
    double *eB;
    int offset_X_hat = 0;
    int offset_A = 0;
    int offset_B = 0;
    int offset_e = 0;
    int offset_U = 0;
    int offset_V = 0;
    int offset_W = 0;
    int offset_eA = 0;
    int offset_eB = 0;
    
    double da[12];
    double db[3];
    
    /* m and n */
    m = mxGetN(pin[0]);
    n = mxGetN(pin[1]);
    
    /* num_a */
    num_a = MAX_NUM_A;
    
    /* X_hat (2xnxm) */
    dim_X_hat[0] = 2; dim_X_hat[1] = n; dim_X_hat[2] = m;
    pout[0] = mxCreateNumericArray(3, dim_X_hat, mxDOUBLE_CLASS, mxREAL);
    X_hat = mxGetPr(pout[0]);
    
    /* A (2 x num_a x n x m) */
    dim_A[0] = 2; dim_A[1] = num_a; dim_A[2] = n; dim_A[3] = m;
    pout[1] = mxCreateNumericArray(4, dim_A, mxDOUBLE_CLASS, mxREAL);
    A = mxGetPr(pout[1]);
    
    /* B (2x3xnxm) */
    dim_B[0] = 2; dim_B[1] = 3; dim_B[2] = n; dim_B[3] = m;
    pout[2] = mxCreateNumericArray(4, dim_B, mxDOUBLE_CLASS, mxREAL);
    B = mxGetPr(pout[2]);
    
    /* e (2xnxm) */
    dim_e[0] = 2; dim_e[1] = n; dim_e[2] = m;
    pout[3] = mxCreateNumericArray(3, dim_e, mxDOUBLE_CLASS, mxREAL);
    e = mxGetPr(pout[3]);
    
    /* U (num_a x num_a x m) */
    dim_U[0] = num_a; dim_U[1] = num_a; dim_U[2] = m;
    pout[4] = mxCreateNumericArray(3, dim_U, mxDOUBLE_CLASS, mxREAL);
    U = mxGetPr(pout[4]);
    
    /* V (3x3xn) */
    dim_V[0] = 3; dim_V[1] = 3; dim_V[2] = n;
    pout[5] = mxCreateNumericArray(3, dim_V, mxDOUBLE_CLASS, mxREAL);
    V = mxGetPr(pout[5]);
    
    /* W (num_a x3xnxm) */
    dim_W[0] = num_a; dim_W[1] = 3; dim_W[2] = n; dim_W[3] = m;
    pout[6] = mxCreateNumericArray(4, dim_W, mxDOUBLE_CLASS, mxREAL);
    W = mxGetPr(pout[6]);
    
    /* eA (num_a x m) */
    pout[7] = mxCreateDoubleMatrix(num_a, m, mxREAL);
    eA = mxGetPr(pout[7]);
    
    /* eB (3xn) */
    pout[8] = mxCreateDoubleMatrix(3, n, mxREAL);
    eB = mxGetPr(pout[8]);
    
    /*-------------------------------------------------------------------
     * Process
     *------------------------------------------------------------------*/

    /* X_hat, A, B, e */
    for ( j = 0 ; j < m ; j ++ )
    {
        for ( i = 0 ; i < n ; i ++ )
        {
            if ( pvisible[i+n*j] )
            {
                /* X_hat */
                reproject_projective_point( pa+num_a*j, pb+3*i, X_hat+offset_X_hat );
                
                /* A */
                for ( k = 0 ; k < num_a ; k ++ )
                {
                    memset(da, 0, num_a*sizeof(double));
                    da[k] = 1;
                    derivative_projective_camera( pa+num_a*j, da, X_hat+offset_X_hat,
                                                  pb+3*i, A+offset_A);
                    offset_A += 2;
                }

                /* B */
                for ( k = 0 ; k < 3 ; k ++ )
                {
                    memset(db, 0, 3*sizeof(double));
                    db[k] = 1;
                    derivative_projective_point( pb+3*i, db, X_hat+offset_X_hat,
                                                 pa+num_a*j, B+offset_B);
                    offset_B += 2;
                }
                
                /* e */
                e[offset_e  ] = pX[offset_X_hat  ] - X_hat[offset_X_hat  ];
                e[offset_e+1] = pX[offset_X_hat+1] - X_hat[offset_X_hat+1];
                offset_e += 2;
            }
            else
            {
                /* X_hat */
                X_hat[offset_X_hat  ] = pX[offset_X_hat  ];
                X_hat[offset_X_hat+1] = pX[offset_X_hat+1];
                
                /* A */
                for ( k = 0 ; k < num_a ; k ++ )
                {
                    A[offset_A  ] = 0;
                    A[offset_A+1] = 0;
                    offset_A += 2;
                }
                
                /* B */
                for ( k = 0 ; k < 3 ; k ++ )
                {
                    B[offset_B  ] = 0;
                    B[offset_B+1] = 0;
                    offset_B += 2;
                }
                
                /* e */
                e[offset_e  ] = 0;
                e[offset_e+1] = 0;
                offset_e += 2;
            }

            offset_X_hat += 2;
        }
    }
        
    /* U, V, W, eA, eB */
    for ( j = 0 ; j < m ; j ++ )
    {
        for ( i = 0 ; i < n ; i ++ )
        {
            /* compute offset */
            offset_U = num_a*num_a*j;
            offset_V = 3*3*i;
            offset_W = num_a*3*(i+n*j);
            offset_A = num_a*2*(i+n*j);
            offset_B = 3*2*(i+n*j);
            offset_eA = num_a*j;
            offset_eB = 3*i;
            offset_e  = 2*(i+n*j);
            
            /* U(:,:,j) = U(:,:,j) + A(:,:,i,j)' * A(:,:,i,j) */
            for ( col = 0 ; col < num_a ; col ++ )
            {
                for ( row = 0 ; row < num_a ; row ++ )
                {
                    U[row+num_a*col+offset_U] += (
                        A[  2*row+offset_A]*A[  2*col+offset_A] +
                        A[1+2*row+offset_A]*A[1+2*col+offset_A]
                    );
                }
            }
            
            /* V(:,:,i) = V(:,:,i) + B(:,:,i,j)' * B(:,:,i,j) */
            for ( col = 0 ; col < 3 ; col ++ )
            {
                for ( row = 0 ; row < 3 ; row ++ )
                {
                    V[row+3*col+offset_V] += (
                        B[  2*row+offset_B]*B[  2*col+offset_B] +
                        B[1+2*row+offset_B]*B[1+2*col+offset_B]
                    );
                }
            }
            
            /* W(:,:,i,j) = A(:,:,i,j)' * B(:,:,i,j) */
            for ( col = 0 ; col < 3 ; col ++ )
            {
                for ( row = 0 ; row < num_a ; row ++ )
                {
                    W[row+num_a*col+offset_W] += (
                        A[  2*row+offset_A]*B[  2*col+offset_B] +
                        A[1+2*row+offset_A]*B[1+2*col+offset_B]
                    );
                }
            }

            /* eA(:,j) = eA(:,j) + A(:,:,i,j)' * e(:,i,j) */
            for ( row = 0 ; row < num_a ; row ++ )
            {
                eA[row+offset_eA] += (
                    A[  2*row+offset_A]*e[  offset_e] +
                    A[1+2*row+offset_A]*e[1+offset_e]
                );
            }

            /* eB(:,i) = eB(:,i) + B(:,:,i,j)' * e(:,i,j) */
            for ( row = 0 ; row < 3 ; row ++ )
            {
                eB[row+offset_eB] += (
                    B[  2*row+offset_B]*e[  offset_e] +
                    B[1+2*row+offset_B]*e[1+offset_e]
                );
            }
        }
    }

    return;
}
