// Amrex
#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_RealBox.H>
#include <AMReX_CoordSys.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Particles.H>
#include <AMReX_Particle.H>
#include <AMReX_Utility.H>
#include <AMReX_ParticleUtil.H>
#include <AMReX_NeighborParticles.H>
#include "AMReX_Array4.H"


// std c++
#include <iostream>
#include <assert.h>     /* assert */
// Own
#include "AMReX_REAL.H"
#include "particle_defs.hpp"
#include "cmath"


void div_F(const amrex::Geometry geom ,amrex::Box const& bx, amrex::Array4<amrex::Real> const& a){
   const auto lo = amrex::lbound(bx);
   const auto hi = amrex::ubound(bx);
   for     (int k = lo.z; k <= hi.z; ++k) {
     for   (int j = lo.y; j <= hi.y; ++j) {
       for (int i = lo.x; i <= hi.x; ++i) { 
           auto res=(a(i+1,j,k,X)-a(i,j,k,X))+
           (a(i,j+1,k,Y)-a(i,j,k,Y))+
           (a(i,j,k+1,Z)-a(i,j,k,Z));

           if(res!=0){
            std::cout << "i,j,k: " << i << ","<< j << "," <<k <<std::endl;
            std::cout << res << std::endl;
           }
       }
     }
   }


}


// Central difference first order


std::array<amrex::Real,3> curl_cdiff_1(amrex::Array4<amrex::Real> const& a ,int i , int j ,int k, const double* ics){
    return {
        ((a(i,j+1,k,Z)-a(i,j-1,k,Z))*ics[Y]*0.5-( a(i,j,k+1,Y)-a(i,j,k-1,Y))*ics[Z]*0.5),
        ((a(i,j,k+1,X)-a(i,j,k-1,X))*ics[Z]*0.5-( a(i+1,j,k,Z)-a(i-1,j,k,Z))*ics[X]*0.5),    
        ((a(i+1,j,k,Y)-a(i-1,j,k,Y))*ics[X]*0.5-( a(i,j+1,k,X)-a(i,j-1,k,X))*ics[Y]*0.5)
    };
}


// forward difference, first order
std::array<amrex::Real,3> curl_fdiff_1(amrex::Array4<amrex::Real> const& a ,int i , int j ,int k, const double* ics){
    return {
        ((a(i,j+1,k,Z)-a(i,j,k,Z))*ics[Y]-( a(i,j,k+1,Y)-a(i,j,k,Y))*ics[Z]),
        ((a(i,j,k+1,X)-a(i,j,k,X))*ics[Z]-( a(i+1,j,k,Z)-a(i,j,k,Z))*ics[X]),
        ((a(i+1,j,k,Y)-a(i,j,k,Y))*ics[X]-( a(i,j+1,k,X)-a(i,j,k,X))*ics[Y])
    };
}

// backward difference, first order
std::array<amrex::Real,3> curl_bdiff_1(amrex::Array4<amrex::Real> const& a ,int i , int j ,int k, const double* ics){
    return {
        ((a(i,j,k,Z)-a(i,j-1,k,Z))*ics[Y]-( a(i,j,k,Y)-a(i,j,k-1,Y))*ics[Z]),
        ((a(i,j,k,X)-a(i,j,k-1,X))*ics[Z]-( a(i,j,k,Z)-a(i-1,j,k,Z))*ics[X]),
        ((a(i,j,k,Y)-a(i-1,j,k,Y))*ics[X]-( a(i,j,k,X)-a(i,j-1,k,X))*ics[Y])
    };
}


void push_B_E(const amrex::Geometry geom, amrex::Box const& bx,  amrex::Array4<amrex::Real> const& B, amrex::Array4<amrex::Real> const& E,double dt){
   const auto lo = amrex::lbound(bx);
   const auto hi = amrex::ubound(bx);
   const auto ics = geom.InvCellSize() ;
   for     (int k = lo.z; k <= hi.z; ++k) {
     for   (int j = lo.y; j <= hi.y; ++j) {
       for (int i = lo.x; i <= hi.x; ++i) { 
        auto curl = curl_cdiff_1(E,i,j,k,ics);
         B(i,j,k,X) -=dt*curl[0];  
         B(i,j,k,Y) -=dt*curl[1];
         B(i,j,k,Z) -=dt*curl[2];
       }
     }
   }

}

void push_E_B(const amrex::Geometry geom, amrex::Box const& bx,  amrex::Array4<amrex::Real> const& E, amrex::Array4<amrex::Real> const& B,double dt){
   const auto lo = amrex::lbound(bx);
   const auto hi = amrex::ubound(bx);
   const auto ics = geom.InvCellSize() ;
   for     (int k = lo.z; k <= hi.z; ++k) {
     for   (int j = lo.y; j <= hi.y; ++j) {
       for (int i = lo.x; i <= hi.x; ++i) { 
        auto curl = curl_cdiff_1(B,i,j,k,ics);
         E(i,j,k,X) +=dt*curl[0];  
         E(i,j,k,Y) +=dt*curl[1];
         E(i,j,k,Z) +=dt*curl[2];
       }
     }
   }

}


void G_Theta_E(const amrex::Geometry geom, amrex::MultiFab&E, amrex::MultiFab&B,double dt ){

    for (amrex::MFIter mfi(E); mfi.isValid(); ++mfi){
        const amrex::Box& box = mfi.validbox();
        amrex::FArrayBox& fab = E[mfi];
        amrex::FArrayBox& fabB = B[mfi];
        amrex::Array4<amrex::Real> const& E_loc = fab.array();
        amrex::Array4<amrex::Real> const& B_loc = fabB.array(); 
        push_B_E(geom,box, B_loc,E_loc,dt);

    }

    B.FillBoundary(geom.periodicity());
}

void G_Theta_B(const amrex::Geometry geom, amrex::MultiFab&E, amrex::MultiFab&B,double dt ){

    for (amrex::MFIter mfi(E); mfi.isValid(); ++mfi){
        const amrex::Box& box = mfi.validbox();
        amrex::FArrayBox& fab = E[mfi];
        amrex::FArrayBox& fabB = B[mfi];
        amrex::Array4<amrex::Real> const& E_loc = fab.array();
        amrex::Array4<amrex::Real> const& B_loc = fabB.array(); 
        push_E_B(geom,box, E_loc,B_loc,dt);

    }

    E.FillBoundary(geom.periodicity());
}



void Theta_2(const amrex::Geometry geom, amrex::MultiFab&E, amrex::MultiFab&B,double dt ){
        G_Theta_E(geom,E,B,dt/2);
        G_Theta_B(geom,E,B,dt);
        G_Theta_E(geom,E,B,dt/2); 
}

template <unsigned int l>
void Theta(const amrex::Geometry geom, amrex::MultiFab&E, amrex::MultiFab&B,double dt ){
    const double alpha=1/(2-std::pow(2,1/(2*l+1)));
    const double beta = 1- 2*alpha;
    assert (l>0);

    Theta<l-1>(geom,E,B,alpha*dt);
    Theta<l-1>(geom,E,B,beta*dt);
    Theta<l-1>(geom,E,B,alpha*dt);
    return;
}

template <>
void Theta<0>(const amrex::Geometry geom, amrex::MultiFab&E, amrex::MultiFab&B,double dt ){
        const double alpha=1/(2-std::pow(2,1/(2*1+1)));
        const double beta = 1- 2*alpha;
        Theta_2(geom,E,B,alpha*dt);
        Theta_2(geom,E,B,beta*dt);
        Theta_2(geom,E,B,alpha*dt);
}



void main_main();

void init_E (const amrex::Geometry geom ,amrex::Box const& bx, amrex::Array4<amrex::Real> const& a)
{
   const auto lo = amrex::lbound(bx);
   const auto hi = amrex::ubound(bx);
   for     (int k = lo.z; k <= hi.z; ++k) {
     for   (int j = lo.y; j <= hi.y; ++j) {
       for (int i = lo.x; i <= hi.x; ++i) { 
         a(i,j,k,X) = 0.0;
         a(i,j,k,Y) = 0.0; 
         a(i,j,k,Z) = 0.0;       
       }
     }
   }
}


void init_B (const amrex::Geometry geom ,amrex::Box const& bx, amrex::Array4<amrex::Real> const& a)
{
   const auto lo = amrex::lbound(bx);
   const auto hi = amrex::ubound(bx);
   for     (int k = lo.z; k <= hi.z; ++k) {
     for   (int j = lo.y; j <= hi.y; ++j) {
       for (int i = lo.x; i <= hi.x; ++i) { 
        if( i % 50 == 0 ){
         a(i,j,k,X) = 0;
         a(i,j,k,Y) = 1;//sin( ( 2*( geom.ProbLo(X) + i*geom.CellSize(X) +1))*3.1415962);
         a(i,j,k,Z) = 0;  
        }
       }
     }
   }

}




int main(int argc, char* argv[])
{

    amrex::Initialize(argc,argv);
    main_main();
    amrex::Finalize();

}
void main_main()
{


    const  int n_cell = 256;
    int max_grid_size=256;
    int nsteps=800;
    double dt=1.0/500;
    // Do a quite even load balancing
    amrex::DistributionMapping::strategy(amrex::DistributionMapping::KNAPSACK);

    // Periodic
    //amrex::Array<int,AMREX_SPACEDIM> is_periodic {AMREX_D_DECL(1,1,1)};
    amrex::Vector<int> is_periodic(AMREX_SPACEDIM,1);     
    // Nodal indexing
    amrex::IndexType typ({AMREX_D_DECL(1,1,1)});


    amrex::IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
    amrex::IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, 8));
    amrex::Box domain(dom_lo, dom_hi,typ);
    amrex::BoxArray ba(domain);

    // Initialize the boxarray "ba" from the single box "bx"
    // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
    ba.maxSize(max_grid_size);

    // This defines the physical box, [-1,1] in each direction.
    amrex::RealBox real_box({AMREX_D_DECL(-1,-1,-1)},
                     {AMREX_D_DECL( 1,1,1)});

    // This defines a Geometry object
    amrex::Geometry geom(domain,&real_box,amrex::CoordSys::cartesian,is_periodic.data());

    // How Boxes are distrubuted among MPI processes
    amrex::DistributionMapping dm(ba);

    int Nghost = 4;
    
    // Ncomp = number of components for each array
    int Ncomp  = 3;
  
    amrex::MultiFab E(ba,dm,Ncomp,Nghost);
    amrex::MultiFab B(ba,dm,Ncomp,Nghost);
    for (amrex::MFIter mfi(E); mfi.isValid(); ++mfi) // Loop over grids
{
    // This is the valid Box of the current FArrayBox.
    // By "valid", we mean the original ungrown Box in BoxArray.
    const amrex::Box& box = mfi.validbox();

    // A reference to the current FArrayBox in this loop iteration.
    amrex::FArrayBox& fab = E[mfi];
    amrex::FArrayBox& fabB = B[mfi];
    // Obtain Array4 from FArrayBox.  We can also do
    //     Array4<Real> const& a = mf.array(mfi);
    amrex::Array4<amrex::Real> const& a = fab.array();
    amrex::Array4<amrex::Real> const& b = fabB.array();

    // Call function f1 to work on the region specified by box.
    // Note that the whole region of the Fab includes ghost
    // cells (if there are any), and is thus larger than or
    // equal to "box".
    init_E(geom,box, a);
    init_B(geom,box, b);


}
E.FillBoundary(geom.periodicity());
B.FillBoundary(geom.periodicity());



for(int step=0; step < nsteps ; step++){
        std::cout << "Step:" << step << std::endl;

    Theta_2(geom,E,B,dt);
    if(step % 5 ==0){
    int n=step;
    amrex::Real time=step*dt;
    const std::string& pltfile_E = amrex::Concatenate("plt_E",n,0);
    WriteSingleLevelPlotfile(pltfile_E, E, {"E_x","E_y","E_z"}, geom, time, n);
    const std::string& pltfile_B = amrex::Concatenate("plt_B",n,0);
    WriteSingleLevelPlotfile(pltfile_B, B, {"B_x","B_y","B_z"}, geom, time, n);
    }
 
}
}