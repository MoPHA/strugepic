#ifndef INTERPOL
#define INTERPOL
#include<array>

namespace strugepic
{ 
    double I_W12(double a,double b);
    double W ( double x , double y , double z);
    double W1 (double x);
    double W12 (double x);
    double W1d (double x);
    double I_W1(double a,double b);
    std::array<double,3> W_s1(const std::array<double,3> coord);

}
#endif
