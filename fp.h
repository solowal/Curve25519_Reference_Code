//
//  fp.h
//  math_ecc
//
//  Created by info on 2016. 9. 4..
//  Copyright (c) 2016년 info. All rights reserved.
//  Developer Hwajeong Seo (https://scholar.google.co.kr/citations?user=NJTYsz0AAAAJ&hl=ko)

#ifndef __math_ecc__fp__
#define __math_ecc__fp__

#include <stdio.h>
#include <iostream>
#include <ctime>
#include <inttypes.h>
#include <stdint.h>

using namespace std;

typedef __uint128_t u128;
typedef uint64_t u64;
typedef uint32_t u32;

int int_add_c99(u32 *r, const u32 *a, const u32 *b, int len);
int int_sub_c99(u32 *r, const u32 *a, const u32 *b, int len);
void int_mul_c99(u32 *r, const u32 *a, const u32 *b, int len);
void int_sqr_c99(u32 *r, const u32 *a, int len);

void fp_add(u32* c, u32* a, u32* b);
void fp_sub(u32* c, u32* a, u32* b);

void mul256_func(u32* c, u32* a, u32* b);
void sqr256_func(u32* c, u32* a);
void red256_func(u32* c, u32* a);

static u32 prod[16]={0,};

#define fp_mul(c,a,b) do{    \
mul256_func(prod,a,b);       \
red256_func(c,prod);         \
}while(0)

#define fp_sqr(c,a) do{      \
sqr256_func(prod,a);         \
red256_func(c,prod);         \
}while(0)

void fp_inv(u32* c, u32* z);

#endif /* defined(__math_ecc__fp__) */
