//
//  ep.h
//  math_ecc
//
//  Created by info on 2016. 9. 14..
//  Copyright (c) 2016년 info. All rights reserved.
//

#ifndef __math_ecc__ep__
#define __math_ecc__ep__

#include "fp.h"
#include <stdio.h>

typedef struct point
{
    uint32_t x[8];
    uint32_t y[8];
    uint32_t z[8];
    uint32_t proj;
} POINT;

void PointAdd(POINT *R, POINT *P, POINT *Q, uint32_t *xd);
void PointDbl(POINT *R, POINT *P);
void RecoverY(POINT* R, POINT* P1, POINT* P2, POINT* P, uint32_t* b);
void ProToAff(POINT* R, POINT* P);
void PointXZMulSecure(POINT* R1, POINT* R2, uint32_t* k, POINT* P);
void PointMul(POINT* R, uint32_t* k, POINT* P, uint32_t* b);

#endif /* defined(__math_ecc__ep__) */
