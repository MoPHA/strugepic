#include"internal_interpolation.hpp"
#ifndef INTERPOL
#define INTERPOL

#define POLY_DEGREE W1_8

// This is just some basic stuff
// Need to think about the usability
namespace strugepic
{    
    // Change this and recompile to have a different interpolating function
    auto const fp_d = &strugepicInternal::POLY_DEGREE<double>;
    
__attribute__ ((visibility ("default")))   double W ( double x , double y , double z){
        return fp_d(x)*fp_d(y)*fp_d(z);
    }

    auto const fp_f = &strugepicInternal::POLY_DEGREE<float>;

__attribute__ ((visibility ("default")))    float W (float x,float y,float z ){
        return fp_f(x)*fp_f(y)*fp_f(z);
    }

__attribute__ ((visibility ("default")))    double W (double* p){
        return fp_f(p[0])*fp_f(p[1])*fp_f(p[2]);
    }

__attribute__ ((visibility ("default")))    float  W (float* p){
        return fp_f(p[0])*fp_f(p[1])*fp_f(p[2]);
    }

}

#endif
