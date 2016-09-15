//
//  fp.c
//  math_ecc
//
//  Created by info on 2016. 9. 4..
//  Copyright (c) 2016년 info. All rights reserved.
//  Developer Hwajeong Seo (https://scholar.google.co.kr/citations?user=NJTYsz0AAAAJ&hl=ko)

#include "fp.h"

// multiprecision addition (32-bit wise)
int int_add_c99(u32 *r, const u32 *a, const u32 *b, int len)
{
    int i;
    u64 sum;
    
    sum = 0;
    for (i = 0; i < len; i++)
    {
        sum += (u64) a[i] + b[i];
        r[i] = (u32) sum;
        sum >>= 32;
    }
    
    return (int) sum; // carry bit
}

// finite field addition (256-bit case)
void fp_add(u32* c, u32* a, u32* b)
{
    int carry;
    u32 mod[] ={0xFFFFFFED,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0x7FFFFFFF};
    
    carry = int_add_c99(c,a,b,8);
    while(carry>0){
        carry -= int_sub_c99(c,c,mod,8);
    }
}

// multiprecision subtraction (32-bit wise)
int int_sub_c99(u32 *r, const u32 *a, const u32 *b, int len)
{
    int i;
    u64 dif;
    
    dif = 1;
    for (i = 0; i < len; i++)
    {
        dif += (u64) a[i] + (~b[i]); // 1's complement http://www.programiz.com/c-programming/bitwise-operators#complement
        r[i] = (u32) dif;
        dif >>= 32;
    }
    
    return (1 - ((int) dif)); // borrow bit
}

// finite field subtraction (256-bit case)
void fp_sub(u32* c, u32* a, u32* b)
{
    int borrow;
    u32 mod[] ={0xFFFFFFED,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0x7FFFFFFF};
    
    borrow = int_sub_c99(c,a,b,8);
    while(borrow>0){
        borrow -= int_add_c99(c,c,mod,8);
    }
}


// multiprecision multiplication (32-bit wise)
void int_mul_c99(u32 *r, const u32 *a, const u32 *b, int len)
{
    int i, j, k;
    u128 accu;
    
    // first iteration of 1st outer loop
    accu = (u64) a[0]*b[0];
    r[0] = (u32) accu;
    accu >>= 32;
    // len-1 iterations of 1st outer loop
    for (i = 1; i <= len-1; i++)
    {
        j = 0; k = i;
        while (k >= 0) accu += (u64) a[j++]*b[k--];
        r[i] = (u32) accu;
        accu >>= 32;
    }
    // len-2 iterations of 2nd outer loop
    for (i = len; i <= 2*len-3; i++)
    {
        k = len-1; j = i-k;
        while (j <= len-1) accu += (u64) a[j++]*b[k--];
        r[i] = (u32) accu;
        accu >>= 32;
    }
    // last iteration of 2nd outer loop
    accu += (u64) a[len-1]*b[len-1];
    r[2*len-2] = (u32) accu;
    r[2*len-1] = (u32) (accu >> 32);
}


// multiprecision squaring (32-bit wise)
void int_sqr_c99(u32 *r, const u32 *a, int len)
{
    int i, j, k;
    u128 accu = 0;
    
    r[0] = 0;
    // len-1 iterations of 1st outer loop
    for(i = 1; i <= len-1; i++)
    {
        j = 0; k = i;
        while (j < k) accu += (u64) a[j++]*a[k--];
        r[i] = (u32) accu;
        accu >>= 32;
    }
    // len-2 iterations of 2nd outer loop
    for(i = len; i <= 2*len-3; i++)
    {
        k = len-1; j = i-k;
        while (j < k) accu += (u64) a[j++]*a[k--];
        r[i] = (u32) accu;
        accu >>= 32;
    }
    r[i] = (u32) accu;
    r[i+1] = 0;
    // double the result obtained so far and add all a[i]*a[i]
    accu = 0;
    for(i = 0; i <= 2*len-2; i += 2)
    {
        k = i>>1;
        accu += (u64) a[k]*a[k] + r[i] + r[i];
        r[i] = (u32) accu;
        accu >>= 32;
        accu += (u64) r[i+1] + r[i+1];
        r[i+1] = (u32) accu;
        accu >>= 32;
    }
}

// multiplication (256-bit case)
void mul256_func(u32* c, u32* a, u32* b){
    int_mul_c99(c, a, b, 8);
}

// squaring (256-bit case)
void sqr256_func(u32* c, u32* a){
    int_sqr_c99(c, a, 8);
}

// reduction (256-bit case)
void red256_func(u32* c, u32* a)
{
    u32 mod[8]={19*2,0,0,0,0,0,0,0};
    u32 tmp1[16]={0,};
    u32 tmp2[16]={0,};
    int carry;
    
    mul256_func(tmp1,&a[8],mod);
    
    carry = int_add_c99(c,a,tmp1,8);
    tmp1[8] = tmp1[8] + carry;
    
    
    mul256_func(tmp2,&tmp1[8],mod);
    
    u32 mod2[] ={0xFFFFFFED,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0x7FFFFFFF};
    
    carry = int_add_c99(c,c,tmp2,8);
    while(carry>0){
        carry -= int_sub_c99(c,c,mod2,8);
    }
    
}

// finite field inversion (256-bit case)
void fp_inv(uint32_t* c, uint32_t* z){
    uint32_t z2[8], z9[8], z11[8], t1[8], t2[8], t3[8], t4[8], t5[8], t6[8];
    int i;
    // z2:=z^(2^1);
    fp_sqr(z2,z);
    
    //z9:=z2^(2^2)*z;
    fp_sqr(z9,z2);
    fp_sqr(z9,z9);
    fp_mul(z9,z9,z);
    
    //z11:=z9*z2;
    fp_mul(z11,z9,z2);
    
    //t1:=z11^(2^1)*z9;
    fp_sqr(t1,z11);
    fp_mul(t1,t1,z9);
    
    //t2:=t1^(2^5)*t1;
    fp_sqr(t2,t1);
    for(i=0;i<4;i++){
        fp_sqr(t2,t2);
    }
    fp_mul(t2,t2,t1);
    
    //t3:=t2^(2^10)*t2;
    fp_sqr(t3,t2);
    for(i=0;i<9;i++){
        fp_sqr(t3,t3);
    }
    fp_mul(t3,t3,t2);
    
    
    //t4:=(t3^(2^20)*t3)^(2^10)*t2;
    fp_sqr(t4,t3);
    for(i=0;i<19;i++){
        fp_sqr(t4,t4);
    }
    fp_mul(t4,t4,t3);
    for(i=0;i<10;i++){
        fp_sqr(t4,t4);
    }
    fp_mul(t4,t4,t2);
    
    //t5:=t4^(2^50)*t4;
    fp_sqr(t5,t4);
    for(i=0;i<49;i++){
        fp_sqr(t5,t5);
    }
    fp_mul(t5,t5,t4);
    
    //t6:=((t5^(2^100)*t5)^(2^50)*t4)^(2^5)*z11;
    fp_sqr(t6,t5);
    for(i=0;i<99;i++){
        fp_sqr(t6,t6);
    }
    fp_mul(t6,t6,t5);
    
    
    for(i=0;i<50;i++){
        fp_sqr(t6,t6);
    }
    fp_mul(t6,t6,t4);
    
    for(i=0;i<5;i++){
        fp_sqr(t6,t6);
    }
    fp_mul(c,t6,z11);
    
}
