#include<math.h>
#include"poly_util.hpp"

#ifndef REAL
#define REAL double
#endif




#ifdef WITH_INTERPOLATION
#ifdef INTERPOLATION_P8R2
extern const int __attribute__((weak)) interpolation_range=2;
// 8th order piecewise polynominal coefficients, from highest to lowest, comments indicate the range
// Function is zero outside the range
constexpr  REAL P8_coeffs[36]={
        15.0/1024,15.0/128,49.0/128,21.0/32,35.0/64,0,0,1,1,  // -2 -> -1
        -15.0/1024,15.0/128,7.0/16,21.0/32,175.0/256,0,-105.0/128,0,337.0/512, // -1 -> 0
        -15.0/1024,-15.0/128,7.0/16,-21.0/32,175.0/256,0,-105.0/128,0,337.0/512, // 0 -> 1
        15.0/1024,-15.0/128,49.0/128,-21.0/32,35.0/64,0,0,-1,1 //1->2
        };


// Derivative for previous
constexpr  REAL P8d_coeffs[32]={
         15.0/128, 105.0/128, 147.0/64,105.0/32,35.0/16,0,0,1, // -2 -> -1 
         -15.0/128,105.0/128,21.0/8,105.0/32,175.0/64,0,-105.0/64,0,  // -1 -> 0
         -15.0/128,-105.0/128,21.0/8,-105.0/32,175.0/64,0,-105.0/64,0,  // 0  -> 1
         15.0/128, -105.0/128, 147.0/64,-105.0/32 ,35.0/16,0,0,-1 // 1 -> 2 
         };


constexpr REAL P8_agregate_coeffs[27]={
        15.0/1024,0,-28.0/1024,0,-70.0/1024,0,420.0/1024,512.0/1024,175.0/1024, //-1->0
        0,120.0/1024,-420.0/1024,672.0/1024,-630.0/1024,0,420.0/1024,512.0/1024,175.0/1024, // 0->1
        -15.0/1024,120.0/1024,-392.0/1024,672.0/1024,-560.0/1024,0,0,1024.0/1024,0 // 1->2                                  
};
constexpr REAL P8d_agregate_coeffs[24]={
        15.0/128,0,-21.0/128,0 ,-35.0/128,0,105.0/128,64.0/128, // -1 -> 0
        0,105.0/128, -315.0/128, 420.0/128, -315.0/128,0, 105.0/128,64.0/128, // 0->1
        -15.0/128, 105.0/128, -147.0/64,+105.0/32 ,-35.0/16,0,0,+1 // 1-2
};


// Integral coefficients (with extra constant so continuous)
constexpr  REAL P8I_coeffs[40]={
        5.0/3072,15.0/1024,7.0/128,7.0/64,7.0/64,0,0,1.0/2,1,7.0/12, // -2 -> 1
        -5./3072,15./1024,1./16,7./64,35./256,0,-35./128,0,337./512,1.0/2 , // -1->0
        -5./3072,-15./1024,1./16,-7./64,35./256,0,-35./128,0,337./512,1.0/2 , // 0->1
        5.0/3072,-15.0/1024,7.0/128,-7.0/64,7.0/64,0,0,-1.0/2,1,5.0/12 // 1 -> 2
        };


REAL __attribute__((weak)) W1(REAL x){
    return P<REAL,8,-2,2>(x,P8_coeffs); 
}

REAL __attribute__((weak)) W(REAL x,REAL y, REAL z){
    return W1(x)*W1(y)*W1(z);
}

REAL __attribute__((weak)) W1d(REAL x){
    return P<REAL,7,-2,2>(x,P8d_coeffs);
}

REAL __attribute__((weak)) I_W1(REAL a,REAL b ){
    return IP<REAL,8,-2,2>(b,P8I_coeffs)-IP<REAL,8,-2,2>(a,P8I_coeffs);    
}

REAL __attribute__((weak)) Wp(REAL x){
    return P<REAL,7,-1,2>(x,P8d_agregate_coeffs);    
}

REAL __attribute__((weak)) I_Wp(REAL a,REAL b){
        return IP<REAL,7,-1,2>(b,P8_agregate_coeffs)-IP<REAL,7,-1,2>(a,P8_agregate_coeffs);
}
#endif


#ifdef INTERPOLATION_PWL

extern const int __attribute__((weak)) interpolation_range=1;

REAL __attribute__((weak)) W1(REAL x){
    const auto fx=fabs(x);
    if(fx >=1 ){
            return 0;
    }
    else{
        return 1-fx;
    }
}

REAL __attribute__((weak)) W(REAL x,REAL y, REAL z){
    return W1(x)*W1(y)*W1(z);
}

REAL __attribute__((weak)) W1d(REAL x){
    if(x <1 && x >0){
        return -1;
    }
    else if( x > -1 &&  x< 0){
        return 1;
    }
    else{
        return 0;
    }

}

REAL __attribute__((weak)) IH_W1(REAL x){
    if(x > 1){
        return 1;
    }
    else if(x < -1){
        return 0;
    }
    else{
        return -(x*fabs(x)-2*x-1)*0.5;
    }
}

REAL __attribute__((weak)) I_W1(REAL a,REAL b ){
    return IH_W1(b)-IH_W1(a); 
}

REAL __attribute__((weak)) Wp(REAL x){
    if(x <1 && x >=0){
        return 1;
    }
    else{
        return 0;
    }
}

REAL IH_wp(REAL x){
    if(x < 0){
    return 0;
    }
    else if(x > 1){
     return 1;       
    }
    else{
    return x;
    }
}

REAL __attribute__((weak)) I_Wp(REAL a,REAL b){
    return IH_wp(b)-IH_wp(a);
}
#endif


#endif
