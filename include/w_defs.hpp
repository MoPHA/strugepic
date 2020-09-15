#ifndef WDEFS
#define WDEFS

#include <AMReX_REAL.H>


// Interface to the interpolation functions used by the alogrithm
// How and where they are implement is up to the user (Just link against a compiled object file or dynamic library )

amrex::Real W1(amrex::Real x);
amrex::Real Wp(amrex::Real x);
amrex::Real I_W1(amrex::Real a,amrex::Real b);
amrex::Real I_Wp(amrex::Real a,amrex::Real b);
extern const int interpolation_range;


#endif
