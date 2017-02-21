/*
 * reproject_point.h
 *
 * Copyright 2008 (C) Taehee Lee
 *
 * This program is part of VLG, available in the terms of the GNU
 * General Public Licenseversion 2.
 */
#ifndef _REPROJECT_POINT_H_
#define _REPROJECT_POINT_H_

#include "vl/rodrigues.h"
#include <stdio.h>
#include <string.h>

void reproject_point(
    const double K[],
    const double a[],
    const double b[],
    const int num_variableK,
    double x[2])
{
    double K_[9];
    double R [9];
    double Rb[3];
    double x_[3];

    /* Construct K */
    memcpy(K_, K, sizeof(double)*9);
    if ( num_variableK == 1 )
    {
        K_[0] = a[6]; /* K_(1,1) */
        K_[4] = a[6]; /* K_(2,2) */
    }
    else if ( num_variableK == 4 )
    {
        K_[0] = a[6]; /* K_(1,1) */
        K_[4] = a[7]; /* K_(2,2) */
        K_[6] = a[8]; /* K_(1,3) */
        K_[7] = a[9]; /* K_(2,3) */
    }
    
    /* Compute the Rotation matrix */
    vl_rodrigues(R, 0, a);
    
    /* x_ = K * ( rodrigues(a(1:3)) * b + a(4:6) ) */
    Rb[0] = R[0]*b[0] + R[3]*b[1] + R[6]*b[2] + a[3];
    Rb[1] = R[1]*b[0] + R[4]*b[1] + R[7]*b[2] + a[4];
    Rb[2] = R[2]*b[0] + R[5]*b[1] + R[8]*b[2] + a[5];
    x_[0] = K_[0]*Rb[0] + K_[3]*Rb[1] + K_[6]*Rb[2];
    x_[1] = K_[1]*Rb[0] + K_[4]*Rb[1] + K_[7]*Rb[2];
    x_[2] = K_[2]*Rb[0] + K_[5]*Rb[1] + K_[8]*Rb[2];

    /* x = x_(1:2) / x_(3) */
    x[0] = x_[0] / x_[2];
    x[1] = x_[1] / x_[2];
}


#endif // _REPROJECT_POINT_H_
