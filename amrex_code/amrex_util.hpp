#include "AMReX_Array.H"
#include "AMReX_Geometry.H"
#include<AMReX.H>
#include <array>
#include <math.h>


// What cell index is a given point in?
// This is equivalent with the index for the "Lower left" corner
std::array<int,3> get_point_cell(const amrex::Geometry geom,const amrex::RealArray& pos){ 
   int idx= floor((pos[0] - geom.ProbLo(0))/geom.CellSize(0));
   int idy= floor((pos[1] - geom.ProbLo(1))/geom.CellSize(1));
   int idz= floor((pos[2] - geom.ProbLo(2))/geom.CellSize(2));
    return {idx,idy,idz};
}
