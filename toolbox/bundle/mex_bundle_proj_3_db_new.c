/*
 * mex_bundle_3_db_new.c
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

void mexFunction( int nout, mxArray *pout[],
                  int nin,  const mxArray *pin[] )
{
    int i, j, row;

    /*
     * get parameters
     * 
     * W:  pin[0]       (num_a x3xnxm)
     * da: pin[1]       (num_a*m x 1)
     * eB: pin[2]       (3xn)
     * Vinv: pin[3]     (3x3xn);
     * a:  pin[4]       (num_a x m)
     * b:  pin[5]       (3xn)
     * X:  pin[6]       (2xnxm)
     * visible: pin[7]  (nxm)
     */
    double *W = mxGetPr(pin[0]);
    double *da = mxGetPr(pin[1]);
    double *eB = mxGetPr(pin[2]);
    double *Vinv = mxGetPr(pin[3]);
    double *a = mxGetPr(pin[4]);
    double *b = mxGetPr(pin[5]);
    double *X = mxGetPr(pin[6]);
    double *visible = mxGetPr(pin[7]);
    int m, n;
    int num_a;

    /*
     * output parameters
     */
    mwSize dim_X_hat[3];
    double *db;
    double *a_new;
    double *b_new;
    double *X_hat;
    int offset_W = 0;
    int offset_da = 0;
    int offset_Vinv = 0;
    int offset_X_hat = 0;
    int offset_db = 0;
    
    double Wda[3];
    
    /* m and n */
    m = mxGetN(pin[4]);
    n = mxGetN(pin[5]);
    
    /* num_a */
    num_a = MAX_NUM_A;

    /* db (3xn) */
    pout[0] = mxCreateDoubleMatrix(3, n, mxREAL);
    db = mxGetPr(pout[0]);
    
    /* a_new (num_a x m) */
    pout[1] = mxCreateDoubleMatrix(num_a, m, mxREAL);
    a_new = mxGetPr(pout[1]);
    
    /* b_new (3xn) */
    pout[2] = mxCreateDoubleMatrix(3, n, mxREAL);
    b_new = mxGetPr(pout[2]);
    
    /* X_hat (2xnxm) */
    dim_X_hat[0] = 2; dim_X_hat[1] = n; dim_X_hat[2] = m;
    pout[3] = mxCreateNumericArray(3, dim_X_hat, mxDOUBLE_CLASS, mxREAL);
    X_hat = mxGetPr(pout[3]);
    
    
    /*-------------------------------------------------------------------
     * Process
     *------------------------------------------------------------------*/

    /* db */
    for ( i = 0 ; i < n ; i ++ )
    {
        /* eB(:,i) */
        Wda[0] = eB[  3*i];
        Wda[1] = eB[1+3*i];
        Wda[2] = eB[2+3*i];
        /* eB(:,i) - sum(W_ij'da_j) */
        for ( j = 0 ; j < m ; j ++ )
        {
            offset_da = num_a*j;
            for ( row = 0 ; row < 3 ; row ++ )
            {
                offset_W = num_a*row + num_a*3*(i+n*j);
                Wda[row] -= (
                    W[  offset_W] * da[  offset_da] +
                    W[1+offset_W] * da[1+offset_da] +
                    W[2+offset_W] * da[2+offset_da] +
                    W[3+offset_W] * da[3+offset_da] +
                    W[4+offset_W] * da[4+offset_da] +
                    W[5+offset_W] * da[5+offset_da]
                );
            }
        }
        /* V_inv(:,:,i) * (eB(:,i) - sum(W_ij'da_j) */
        offset_db = 3*i;
        offset_Vinv = 9*i;
        for ( row = 0 ; row < 3 ; row ++ )
        {
            db[row+offset_db] = (
                Vinv[row  +offset_Vinv] * Wda[0] +
                Vinv[row+3+offset_Vinv] * Wda[1] +
                Vinv[row+6+offset_Vinv] * Wda[2]
            );
        }
    }
    
    /* a_new */
    for ( j = 0 ; j < num_a*m ; j ++ )
    {
        a_new[j] = a[j] + da[j];
    }
    
    /* b_new */
    for ( j = 0 ; j < 3*n ; j ++ )
    {
        b_new[j] = b[j] + db[j];
    }
    
    /* X_hat */
    for ( j = 0 ; j < m ; j ++ )
    {
        for ( i = 0 ; i < n ; i ++ )
        {
            if ( visible[i+n*j] )
            {
                /* X_hat */
                reproject_projective_point( a_new+num_a*j, b_new+3*i, X_hat+offset_X_hat );
            }
            else
            {
                /* X_hat */
                X_hat[offset_X_hat  ] = X[offset_X_hat  ];
                X_hat[offset_X_hat+1] = X[offset_X_hat+1];
            }
            offset_X_hat += 2;
        }
    }
    

    return;
}
