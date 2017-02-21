/*
 * mex_bundle_2_Se.c
 *
 * Copyright 2008 (C) Taehee Lee
 *
 * This program is part of VLG, available in the terms of the GNU
 * General Public Licenseversion 2.
 */
#include <mexutils.h>
#include <stdio.h>
#include <string.h>

#define MAX_NUM_A   10

void mexFunction( int nout, mxArray *pout[],
                  int nin,  const mxArray *pin[] )
{
    int i, j, k, row, col;

    /*
     * get parameters
     * Y: pin[0]        (num_a x3xnxm)
     * W: pin[1]        (num_a x3xnxm)
     * U: pin[2]        (num_a x num_a x m)
     * eA: pin[3]       (num_a xm)
     * eB: pin[4]       (3xn)
     */
    double *Y = mxGetPr(pin[0]);
    double *W = mxGetPr(pin[1]);
    double *U = mxGetPr(pin[2]);
    double *eA = mxGetPr(pin[3]);
    double *eB = mxGetPr(pin[4]);
    int m, n;
    int num_a;

    /*
     * output parameters
     */
    double *S;
    double *e;
    int offset_S = 0;
    int offset_e = 0;
    int offset_Y = 0;
    int offset_U = 0;
    int offset_W = 0;
    int offset_eA = 0;
    int offset_eB = 0;
    
    double Sjk[MAX_NUM_A][MAX_NUM_A];
    double YeB[MAX_NUM_A];
    
    /* m and n */
    m = mxGetN(pin[3]);
    n = mxGetN(pin[4]);
    
    /* num_a */
    num_a = mxGetM(pin[0]);

    /* S (num_a*m x num_a*m) */
    pout[0] = mxCreateDoubleMatrix(num_a*m, num_a*m, mxREAL);
    S = mxGetPr(pout[0]);
    
    /* e (num_a*m x 1) */
    pout[1] = mxCreateDoubleMatrix(num_a*m, 1, mxREAL);
    e = mxGetPr(pout[1]);
    
    /*-------------------------------------------------------------------
     * Process
     *------------------------------------------------------------------*/

    /* S */
    for ( k = 0 ; k < m ; k ++ )
    {
        for ( j = 0 ; j < m ; j ++ )
        {
            /* compute offset */
            offset_S = num_a*(j + num_a*m*k);
            
            /* S_jj = U_(j) */
            if ( j == k )
            {
                offset_U = num_a*num_a*j;
                for ( col = 0 ; col < num_a ; col ++ )
                {
                    for ( row = 0 ; row < num_a ; row ++ )
                    {
                        Sjk[row][col] = U[row+num_a*col+offset_U];
                    }
                }
            }
            else
            {
                offset_U = num_a*num_a*j;
                for ( col = 0 ; col < num_a ; col ++ )
                {
                    for ( row = 0 ; row < num_a ; row ++ )
                    {
                        Sjk[row][col] = 0;
                    }
                }
            }
            /* S_jk = S_jk - sum_i( Y_ij W_ik' ) */
            for ( i = 0 ; i < n ; i ++ )
            {
                offset_Y = num_a*3*(i+n*j);
                offset_W = num_a*3*(i+n*k);
                for ( col = 0 ; col < num_a ; col ++ )
                {
                    for ( row = 0 ; row < num_a ; row ++ )
                    {
                        Sjk[row][col] -= (
                            Y[row        +offset_Y]*W[col        +offset_W] +
                            Y[row+num_a  +offset_Y]*W[col+num_a  +offset_W] +
                            Y[row+num_a*2+offset_Y]*W[col+num_a*2+offset_W]
                        );
                    }
                }
            }
            
            /* assign to S */
            for ( col = 0 ; col < num_a ; col ++ )
            {
                for ( row = 0 ; row < num_a ; row ++ )
                {
                    S[row+num_a*m*col+offset_S] = Sjk[row][col];
                }
            }
        }
    }
    
    /* e_ */
    for ( j = 0 ; j < m ; j ++ )
    {
        offset_e  = num_a*j;
        offset_eA = num_a*j;
        memset(YeB, 0, sizeof(double)*num_a);
        for ( i = 0 ; i < n ; i ++ )
        {
            offset_Y  = num_a*3*(i+n*j);
            offset_eB = 3*i;
            for ( row = 0 ; row < num_a ; row ++ )
            {
                YeB[row] += (
                    Y[row        +offset_Y] * eB[  offset_eB] +
                    Y[row+num_a  +offset_Y] * eB[1+offset_eB] +
                    Y[row+num_a*2+offset_Y] * eB[2+offset_eB]
                );
            }
        }
        
        for ( row = 0 ; row < num_a ; row ++ )
        {
            e[row+offset_e] = eA[row+offset_eA] - YeB[row];
        }
    }

    return;
}
