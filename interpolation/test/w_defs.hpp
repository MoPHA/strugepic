#ifndef WDEFS
#define WDEFS



// Interface to the interpolation functions used by the alogrithm
// How and where they are implement is up to the user (Just link against a compiled object file or dynamic library )

double W(double x,double y ,double z);
double W1d(double x);
double W1(double x);
double Wp(double x);
double I_W1(double a,double b);
double I_Wp(double a,double b);


#endif
