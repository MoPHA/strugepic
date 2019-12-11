#include "AMReX_Array.H"
#include "AMReX_Geometry.H"
#include "AMReX_REAL.H"
#include <math.h>
#include "particle_defs.hpp"


// What cell index is a given point in?
// This is equivalent with the index for the "Lower left" corner
amrex::IntArray get_point_cell(const amrex::Geometry geom,const amrex::RealArray pos){ 
   int idx= floor((pos[0] - geom.ProbLo(0))/geom.CellSize(0));
   int idy= floor((pos[1] - geom.ProbLo(1))/geom.CellSize(1));
   int idz= floor((pos[2] - geom.ProbLo(2))/geom.CellSize(2));
    return {idx,idy,idz};
}

// Create a single particle in the simulation (about in the middle)
// This is just used for debugging and some cyclotron stuff

 void add_single_particle( CParticleTile&particlet ,amrex::RealArray pos , amrex::RealArray vel, double m,double q){ 
    CParticle p;
    p.id()   = CParticle::NextID();
    p.cpu()  = amrex::ParallelDescriptor::MyProc();
    p.pos(0) = pos[0];
    p.pos(1) = pos[1];
    p.pos(2) = pos[2];
    p.rdata(0)=m;
    p.rdata(1)=q;
    p.rdata(2)=vel[0];
    p.rdata(3)=vel[1];
    p.rdata(4)=vel[2];
    particlet.push_back(p);
}
