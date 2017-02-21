/*
 * mex_bundle_3_db_new.c
 *
 * Copyright 2008 (C) Taehee Lee
 *
 * This program is part of VLG, available in the terms of the GNU
 * General Public Licenseversion 2.
 */
#include <mexutils.h>
#include "reproject_point.h"

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
     * K:  pin[4]       (3x3)
     * a:  pin[5]       (num_a x m)
     * b:  pin[6]       (3xn)
     * X:  pin[7]       (2xnxm)
     * visible: pin[8]  (nxm)
     */
    double *W = mxGetPr(pin[0]);
    double *da = mxGetPr(pin[1]);
    double *eB = mxGetPr(pin[2]);
    double *Vinv = mxGetPr(pin[3]);
    double *pK = mxGetPr(pin[4]);
    double *a = mxGetPr(pin[5]);
    double *b = mxGetPr(pin[6]);
    double *X = mxGetPr(pin[7]);
    double *visible = mxGetPr(pin[8]);
    double **K;
    int m, n;
    int num_a, num_variableK;

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
    m = mxGetN(pin[5]);
    n = mxGetN(pin[6]);
    
    /* num_a and num_variableK */
    num_a = mxGetM(pin[0]);
    num_variableK = num_a - 6;

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

    /* K */
    K = (double **)malloc(sizeof(double *) * m);
    for ( j = 0 ; j < m ; j ++ )
    {
        K[j] = (double *)malloc(sizeof(double) * 9);
        K[j][0] = pK[0+j*4];   K[j][3] = 0;           K[j][6] = pK[2+j*4];
        K[j][1] = 0;           K[j][4] = pK[1+j*4];   K[j][7] = pK[3+j*4];
        K[j][2] = 0;           K[j][5] = 0;           K[j][8] = 1;
    }
    
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
                reproject_point( K[j], a_new+num_a*j, b_new+3*i, num_variableK, X_hat+offset_X_hat );
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
