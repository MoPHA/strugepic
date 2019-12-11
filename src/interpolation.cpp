#include"internal_interpolation.hpp"
#include<array>

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


__attribute__ ((visibility ("default"))) double W12(double  x){
                    return strugepicInternal::W12<double,POLY_DEGREE>(x);
                }
__attribute__ ((visibility ("default"))) double W1(double x){
                    return strugepicInternal::W1<double,POLY_DEGREE>(x);
}

__attribute__ ((visibility ("default"))) double W1d(double  x){
                    return strugepicInternal::W1d<double,POLY_DEGREE>(x);
                }

 __attribute__ ((visibility ("default")))double I_W12(double a,double b){
                return strugepicInternal::I_W12<double,POLY_DEGREE>(a,b);
                }

 __attribute__ ((visibility ("default")))double I_W1(double a,double b){
                return strugepicInternal::I_W1<double,POLY_DEGREE>(a,b);
                }
 __attribute__ ((visibility ("default"))) std::array<double,3> W_s1(const std::array<double,3> coord){
                return {strugepicInternal::W12<double,3>(coord[0])*strugepicInternal::W1<double,3>(coord[1])*strugepicInternal::W1<double,3>(coord[2]),   
                        strugepicInternal::W1<double,3>(coord[0])*strugepicInternal::W12<double,3>(coord[1])*strugepicInternal::W1<double,3>(coord[2]),
                        strugepicInternal::W1<double,3>(coord[0])*strugepicInternal::W1<double,3>(coord[1])*strugepicInternal::W12<double,3>(coord[2])
                        };
                }

}

