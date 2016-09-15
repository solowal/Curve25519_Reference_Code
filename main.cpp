//
//  main.cpp
//  math_ecc
//
//  Created by info on 2016. 8. 31..
//  Copyright (c) 2016년 info. All rights reserved.
//  Developer Hwajeong Seo (https://scholar.google.co.kr/citations?user=NJTYsz0AAAAJ&hl=ko)


#include "ep.h"

u32 base_x[] = {0xABD6E649,0x89EF9E76,0x6ECA9748,0x9F4B6709,0xD3177BDD,0x4A140E05,0x2ACA6E15,0x445B88B0};
u32 base_y[] = {0x3607448E,0xBD71B8CF,0xB82DF54F,0xD49C10F0,0xDF43B91C,0x3058B89F,0xFC8214D5,0x3ED3CFAA};
u32 base_z[] = {0x1,0x0,0x0,0x0,0x0,0x0,0x0,0x0};
u32 k_in[]   = {0xB74B39CF,0xD4D0042,0x82F3DAEA,0xCAA46734,0x5DBDD571,0xC649CEF6,0x81896E34,0x594B0DCB};//test
u32 b_in[]   = {0xFFF892E5,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0x7FFFFFFF};

int main(int argc, const char * argv[]) {
    POINT P, R;
    
    memcpy(P.x,base_x,sizeof(base_x));
    memcpy(P.y,base_y,sizeof(base_y));
    memcpy(P.z,base_z,sizeof(base_z));
    
    PointMul(&R, k_in, &P, b_in);   //scalar multiplication test
    
    printf("X: \n");
    for (int i=0; i<8; i++) {
        printf("0x%x\n",R.x[i]);
    }

    printf("Y: \n");
    for (int i=0; i<8; i++) {
        printf("0x%x\n",R.y[i]);
    }
    
    printf("Z: \n");
    for (int i=0; i<8; i++) {
        printf("0x%x\n",R.z[i]);
    }
    
    return 0;
}
