#include "AMReX_Array.H"
#include "AMReX_Geometry.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_ParallelReduce.H"
#include "AMReX_REAL.H"
#include <iostream>
#include <math.h>
#include <utility>
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
            std::cout << "(" << p.cpu()<<"," << p.id()<<")" << std::endl;
            std::cout << "POS: [" << p.pos(0)<<"," << p.pos(1)<<"," << p.pos(2) << "]" << std::endl;
            std::cout << "VEL: [" << p.rdata(2)<<"," << p.rdata(3)<<"," << p.rdata(4) << "]" << std::endl;
        }
    }

}

std::pair<amrex::Real,amrex::Real> get_total_energy(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab &E, amrex::MultiFab &B ){
    // Kinetic part 
    amrex::Real E_kin=0;
    amrex::Real E_field=0;
    for (CParIter pti(P, 0); pti.isValid(); ++pti) {
        auto&  particles = pti.GetArrayOfStructs();
        for(auto p : particles ){
            E_kin+=p.rdata(M)*0.5*( p.rdata(VX)*p.rdata(VX)+p.rdata(VY)*p.rdata(VY)+p.rdata(VZ)*p.rdata(VZ)  );
        }
    }
    for (amrex::MFIter mfi(E); mfi.isValid(); ++mfi){
        const amrex::Box& box = mfi.validbox();
        amrex::FArrayBox& fab = E[mfi];
        amrex::FArrayBox& fabB = B[mfi];
        amrex::Array4<amrex::Real> const& E_loc = fab.array();
        amrex::Array4<amrex::Real> const& B_loc = fabB.array(); 

        const auto lo = amrex::lbound(box);
        const auto hi = amrex::ubound(box);
        for     (int k = lo.z; k <= hi.z; ++k) {
            for   (int j = lo.y; j <= hi.y; ++j) {
                for (int i = lo.x; i <= hi.x; ++i) { 
                  E_field+=0.5*( E_loc(i,j,k,X)*E_loc(i,j,k,X)+E_loc(i,j,k,Y)*E_loc(i,j,k,Y)+E_loc(i,j,k,Z)*E_loc(i,j,k,Z));
                  E_field+=0.5*( B_loc(i,j,k,X)*B_loc(i,j,k,X)+B_loc(i,j,k,Y)*B_loc(i,j,k,Y)+B_loc(i,j,k,Z)*B_loc(i,j,k,Z));
                }
            }
        }

    }

    amrex::ParallelAllReduce::Sum(E_kin,amrex::ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Sum(E_field,amrex::ParallelDescriptor::Communicator());
    return std::make_pair(E_field*geom.CellSize(X)*geom.CellSize(Y)*geom.CellSize(Z),E_kin);

}



//
// From
// https://stackoverflow.com/questions/4633177/c-how-to-wrap-a-float-to-the-interval-pi-pi

double wrapMax(double x, double max)
{
    /* integer math: `(max + x % max) % max` */
    return fmod(max + fmod(x, max), max);
}
/* wrap x -> [min,max) */
double wrapMinMax(double x, double min, double max)
{
    return min + wrapMax(x - min, max - min);
}
