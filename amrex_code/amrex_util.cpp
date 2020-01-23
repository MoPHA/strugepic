#include "AMReX_Array.H"
#include "AMReX_CudaContainers.H"
#include "AMReX_Geometry.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_ParallelReduce.H"
#include "AMReX_REAL.H"
#include <iostream>
#include <math.h>
#include <utility>
#include <vector>
#include "particle_defs.hpp"


// Init external dirichle boundary condition

// Either periodic or dirichle

void FillDirichletBoundary(const amrex::Geometry geom, amrex::MultiFab &A,const amrex::Real b_val){
   
   const auto domain=geom.Domain();
   const auto lo = amrex::lbound(domain);
   const auto hi = amrex::ubound(domain);
   const auto periodic=geom.isPeriodic();
    for (amrex::MFIter mfi(A); mfi.isValid(); ++mfi){
        const amrex::Box& box = mfi.fabbox();
        amrex::Array4<amrex::Real> const& a = A.array(mfi); 

    amrex::ParallelFor(box,  [=] AMREX_GPU_DEVICE (int i,int j,int k ){
            if(
              (i < lo.x && !periodic[X] ) ||
              (i > hi.x && !periodic[X] ) ||
              (j < lo.y && !periodic[Y] ) ||
              (j > hi.y && !periodic[Y] ) ||
              (k < lo.z && !periodic[Z] ) ||
              (k > hi.z && !periodic[Z] )  
              ){

                a(i,j,k,X)=b_val;
                a(i,j,k,Y)=b_val;
                a(i,j,k,Z)=b_val;

            }
         });

    }
}




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





// Create a single particle in the simulation

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


int x2_dist(const amrex::Geometry geom,std::array<double,3> coord){

    return (int) -(coord[X]-geom.ProbLo(X))*(coord[X]-geom.ProbHi(X))*50;

}

void add_particle_density(const amrex::Geometry geom , CParticleContainer&P, int (*dist_func)(const amrex::Geometry,double,double,double) ,double m, double q ){

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> distx(0,geom.CellSize(X));
    std::uniform_real_distribution<double> disty(0,geom.CellSize(Y));
    std::uniform_real_distribution<double> distz(0,geom.CellSize(Z));


for(amrex::MFIter mfi= P.MakeMFIter(0) ;mfi.isValid();++mfi){
    
    // Each grid,tile has a their own local particle container
    auto& particles = P.GetParticles(0)[std::make_pair(mfi.index(),
                                        mfi.LocalTileIndex())];
    auto box=mfi.validbox();


   const auto lo = amrex::lbound(box);
   const auto hi = amrex::ubound(box);
   for     (int k = lo.z; k <= hi.z; ++k) {
     for   (int j = lo.y; j <= hi.y; ++j) {
       for (int i = lo.x; i <= hi.x; ++i) { 
           double x = geom.ProbLo(X) + i*geom.CellSize(X);
           double y = geom.ProbLo(Y) + j*geom.CellSize(Y);
           double z = geom.ProbLo(Z) + k*geom.CellSize(Z);
            add_single_particle(particles,{x,y,z},{distx(mt),disty(mt),distz(mt) },m,q);
       }
     }
   }
    }

}


void add_particle_one_per_cell(const amrex::Geometry geom, CParticleContainer&P,double m,double q){
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(-0.01,0.01);


for(amrex::MFIter mfi= P.MakeMFIter(0) ;mfi.isValid();++mfi){
    
    // Each grid,tile has a their own local particle container
    auto& particles = P.GetParticles(0)[std::make_pair(mfi.index(),
                                        mfi.LocalTileIndex())];
    auto box=mfi.validbox();


   const auto lo = amrex::lbound(box);
   const auto hi = amrex::ubound(box);
   for     (int k = lo.z; k <= hi.z; ++k) {
     for   (int j = lo.y; j <= hi.y; ++j) {
       for (int i = lo.x; i <= hi.x; ++i) { 
           double x = geom.ProbLo(X) + (i+0.5)*geom.CellSize(X);
           double y = geom.ProbLo(Y) + (j+0.5)*geom.CellSize(Y);
           double z = geom.ProbLo(Z) + (k+0.5)*geom.CellSize(Z);
            add_single_particle(particles,{x,y,z},{dist(mt),dist(mt),dist(mt) },m,q);
       }
     }
   }
    }
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

    auto E_L2_norm=E.norm2({X,Y,Z});
    auto B_L2_norm=B.norm2({X,Y,Z});
    E_field+=E_L2_norm[X]*E_L2_norm[X]+E_L2_norm[Y]*E_L2_norm[Y]+E_L2_norm[Z]*E_L2_norm[Z];
    E_field+=B_L2_norm[X]*B_L2_norm[X]+B_L2_norm[Y]*B_L2_norm[Y]+B_L2_norm[Z]*B_L2_norm[Z];
    E_field*=0.5;

    amrex::ParallelAllReduce::Sum(E_kin,amrex::ParallelDescriptor::Communicator());
    return std::make_pair(E_field*geom.CellSize(X)*geom.CellSize(Y)*geom.CellSize(Z),E_kin);

}

std::pair<std::array<amrex::Real,3>,std::array<amrex::Real,3>> get_total_momentum(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab &E, amrex::MultiFab &B){
    
    std::array<amrex::Real,3> P_part={0,0,0};
    std::array<amrex::Real,3> P_field={0,0,0};
    for (CParIter pti(P, 0); pti.isValid(); ++pti) {
        auto&  particles = pti.GetArrayOfStructs();
        for(auto p : particles ){
            P_part[X]+=p.rdata(M)*p.rdata(VX);
            P_part[Y]+=p.rdata(M)*p.rdata(VY);
            P_part[Z]+=p.rdata(M)*p.rdata(VZ);
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
                P_field[X]+=E_loc(i,j,k,Y)*B_loc(i,j,k,Z)-E_loc(i,j,k,Z)*B_loc(i,j,k,Y);
                P_field[Y]+=E_loc(i,j,k,Z)*B_loc(i,j,k,X)-E_loc(i,j,k,X)*B_loc(i,j,k,Z);
                P_field[Z]+=E_loc(i,j,k,X)*B_loc(i,j,k,Y)-E_loc(i,j,k,Y)*B_loc(i,j,k,X);
                }
            }
        }

    }
    P_field[X]*=geom.CellSize(X)*geom.CellSize(Y)*geom.CellSize(Z);
    P_field[Y]*=geom.CellSize(X)*geom.CellSize(Y)*geom.CellSize(Z);
    P_field[Z]*=geom.CellSize(X)*geom.CellSize(Y)*geom.CellSize(Z);
    return std::make_pair(P_field,P_part);

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
