#include"internal_interpolation.hpp"
#ifndef INTERPOL
#define INTERPOL

#define POLY_DEGREE 8

// This is just some basic stuff
// Need to think about the usability
namespace strugepic
{    
    // Change this and recompile to have a different interpolating function
    auto const fp_d = &strugepicInternal::W1<double,POLY_DEGREE>;
    
__attribute__ ((visibility ("default")))   double W ( double x , double y , double z){
        return fp_d(x)*fp_d(y)*fp_d(z);
    }

    auto const fp_f = &strugepicInternal::W1<float,POLY_DEGREE>;

__attribute__ ((visibility ("default")))    float W (float x,float y,float z ){
        return fp_f(x)*fp_f(y)*fp_f(z);
    }

__attribute__ ((visibility ("default"))) double W12(double  x){
                    return strugepicInternal::W12<double,POLY_DEGREE>(x);
                }

__attribute__ ((visibility ("default"))) float W12(float  x){
                    return strugepicInternal::W12<float,POLY_DEGREE>(x);
                }

__attribute__ ((visibility ("default"))) double W1d(double  x){
                    return strugepicInternal::W1d<double,POLY_DEGREE>(x);
                }

__attribute__ ((visibility ("default"))) float W1d(float  x){
                    return strugepicInternal::W1d<float,POLY_DEGREE>(x);
                }
}

#endif
