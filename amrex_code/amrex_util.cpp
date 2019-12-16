#include "AMReX_Array.H"
#include "AMReX_Geometry.H"
#include "AMReX_REAL.H"
#include <math.h>
#include <vector>
#include "particle_defs.hpp"


// What cell index is a given point in?
// This is equivalent with the index for the "Lower left" corner
// Vectorization? 
std::array<int,3> get_point_cell(const amrex::Geometry geom,const amrex::RealArray pos){ 
    std::array<int,3> id;
    auto icellsize=geom.InvCellSize();
    auto problo=geom.ProbLo();
    id[0]=floor((pos[0] -problo[0])*icellsize[0]);
    id[1]=floor((pos[1] -problo[1])*icellsize[1]);
    id[2]=floor((pos[2] -problo[2])*icellsize[2]);
    return id;
}


// How many line segments are there if the path is split at each cell boundary
// Exact border case is not handled?
// I.e x_start or x_end at cell boundary
// Vectorization? 
std::array<int,3> get_num_segments(const amrex::Geometry geom,const amrex::RealArray x_start,const amrex::RealArray  x_end){
    std::array<int,3> n;
    auto start_idx=get_point_cell(geom,x_start);
    auto end_idx=get_point_cell(geom,x_end);
    n[0]=abs(start_idx[0]-end_idx[0])+1;
    n[1]=abs(start_idx[1]-end_idx[1])+1;
    n[2]=abs(start_idx[2]-end_idx[2])+1;
    return n;
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

void print_Particle_info(const amrex::Geometry geom,CParticleContainer&P ){

    for (CParIter pti(P, 0); pti.isValid(); ++pti) {
        auto&  particles = pti.GetArrayOfStructs();
        for(auto p : particles ){
            amrex::Print() << "(" << p.cpu()<<"," << p.id()<<")" << std::endl;
            amrex::Print() << "[" << p.pos(0)<<"," << p.pos(1)<<"," << p.pos(2) << "]" << std::endl;
            amrex::Print() << "[" << p.rdata(2)<<"," << p.rdata(3)<<"," << p.rdata(4) << "]" << std::endl;
        }
    }

}
