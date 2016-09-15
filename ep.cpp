//
//  ep.c
//  math_ecc
//
//  Created by info on 2016. 9. 14..
//  Copyright (c) 2016년 info. All rights reserved.
//  Developer Hwajeong Seo (https://scholar.google.co.kr/citations?user=NJTYsz0AAAAJ&hl=ko)

#include "ep.h"


#define get_bit(k, i) (( k[ (i) >> 5 ] >> ((i) & 0x1F)) & 1)

// point addition
void PointAdd(POINT *R, POINT *P, POINT *Q, u32 *xd){
    u32* xp = P->x;
    u32* zp = P->z;
    
    u32* xq = Q->x;
    u32* zq = Q->z;
    
    u32* xr = R->x;
    u32* zr = R->z;
    
    u32 t1[8];
    u32 t2[8];
    
    fp_add(t1,xp,zp);   //t1 := xp+zp;
    fp_sub(t2,xp,zp);   //t2 := xp-zp;
    fp_sub(xr,xq,zq);   //xr := xq-zq;
    fp_mul(zr,t1,xr);   //zr := t1*xr;
    fp_add(t1,xq,zq);   //t1 := xq+zq;
    fp_mul(xr,t1,t2);   //xr := t1*t2;
    fp_sub(t1,xr,zr);   //t1 := xr-zr;
    fp_add(t2,xr,zr);   //t2 := xr+zr;
    fp_sqr(xr,t2);      //xr := t2*t2;
    fp_sqr(t2,t1);      //t2 := t1*t1;
    fp_mul(zr,xd,t2);   //zr := xd*t2;
}

// point doubling
void PointDbl(POINT *R, POINT *P){
    u32* xp = P->x;
    u32* zp = P->z;
    
    u32* xr = R->x;
    u32* zr = R->z;
    
    u32 t1[8];
    u32 t2[8];
    
    u32 c[8] = {0x1DB42, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0, 0x0};
    
    fp_add(t1,xp,zp);   //t1 := xp+zp;
    fp_sqr(t2,t1);      //t2 := t1*t1;
    fp_sub(t1,xp,zp);   //t1 := xp-zp;
    fp_sqr(zr,t1);      //zr := t1*t1;
    fp_mul(xr,t2,zr);   //xr := t2*zr;
    fp_sub(t1,t2,zr);   //t1 := t2-zr;

    fp_mul(t2,t1,c);    //t2 := t1*c;
    
    
    fp_add(t2,t2,zr);   //t2 := t2+zr;
    fp_mul(zr,t1,t2);   //zr := t1*t2;
}

// point recovery
void RecoverY(POINT* R, POINT* P1, POINT* P2, POINT* P, u32* b){
    u32* x1 = P1->x;
    u32* z1 = P1->z;
    
    u32* x2 = P2->x;
    u32* z2 = P2->z;
    
    u32* x = P->x;
    u32* y = P->y;
    
    u32* xr = R->x;
    u32* yr = R->y;
    u32* zr = R->z;
    
    u32 t1[8];
    u32 t2[8];
    u32 t3[8];
    u32 t4[8];
    fp_mul(t1,x,x1);    //t1 := x*x1;
    fp_sub(t1,t1,z1);   //t1 := t1-z1;
    fp_mul(t2,z1,x);    //t2 := z1*x;
    fp_sub(t2,x1,t2);   //t2 := x1-t2;
    fp_mul(t3,z2,t1);   //t3 := z2*t1;
    fp_mul(t4,x2,t2);   //t4 := x2*t2;
    fp_add(t2,b,b);     //t2 := 4*b;
    fp_add(t2,t2,t2);
    fp_mul(t2,t2,y);    //t2 := t2*y;
    fp_mul(t2,t2,z2);   //t2 := t2*z2;
    fp_mul(t2,t2,x2);   //t2 := t2*x2;
    fp_mul(t2,t2,z1);   //t2 := t2*z1;
    fp_mul(zr,t2,z1);   //zr := t2*z1;
    fp_mul(xr,t2,x1);   //xr := t2*x1;
    fp_add(t2,t3,t4);   //t2 := t3+t4;
    fp_sub(t3,t3,t4);   //t3 := t3-t4;
    fp_mul(yr,t2,t3);   //yr := t2*t3;
}

// point from projective to affine
void ProToAff(POINT* R, POINT* P){
    u32 zero[8]={0,};
    u32* xp = P->x;
    u32* yp = P->y;
    u32* zp = P->z;
    
    u32* xr = R->x;
    u32* yr = R->y;
    u32* zr = R->z;
    
    u32 t1[8];
    
    fp_inv(t1,zp);      //t1 := 1/zp;
    fp_mul(yr,yp,t1);   //yr := yp*t1;
    fp_mul(xr,xp,t1);   //xr := xp*t1;
    
    fp_add(xr,xr,zero);
    fp_add(yr,yr,zero);
    memcpy(zr, zero, sizeof(zero));
    zr[0] = 1;
}

// montgomery ladder
void PointXZMulSecure(POINT* R1, POINT* R2, u32* k, POINT* P){
    u32 xp[8];
    memcpy(xp, P->x, sizeof(xp));
    
    POINT   T[2];
    
    memcpy(T[0].x, xp, sizeof(xp));
    memset(T[0].y, 0, sizeof(xp));
    memset(T[0].z, 0, sizeof(xp));
    T[0].z[0] = 1;
    
    PointDbl(&T[1], &T[0]);
    
    for (int i=253; i>-1; i--) {
        int ki = get_bit(k,i);
        PointAdd(&T[1-ki], &T[1-ki], &T[ki], xp);
        PointDbl(&T[ki], &T[ki]);
    }
    memcpy(R1, &T[0], sizeof(POINT));
    memcpy(R2, &T[1], sizeof(POINT));
}

// scalar multiplication
void PointMul(POINT* R, u32* k, POINT* P, u32* b){
    POINT P1, P2;
    PointXZMulSecure(&P1, &P2, k, P);
    RecoverY(&P1, &P1, &P2, P, b);
    ProToAff(R, &P1);
}