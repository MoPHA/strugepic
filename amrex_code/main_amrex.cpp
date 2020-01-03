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
// Own
#include "amrex_util.hpp"
#include "propagators.hpp"
#include "particle_defs.hpp"
#include "cmath"

void main_main();



void init_E (const amrex::Geometry geom ,amrex::Box const& bx, amrex::Array4<amrex::Real> const& a)
{
   const auto lo = amrex::lbound(bx);
   const auto hi = amrex::ubound(bx);
   for     (int k = lo.z; k <= hi.z; ++k) {
     for   (int j = lo.y; j <= hi.y; ++j) {
       for (int i = lo.x; i <= hi.x; ++i) { 
         a(i,j,k,0) = 0;
         a(i,j,k,1) = 0; 
         a(i,j,k,2) = 0;
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
         a(i,j,k,0) = 0.0;
         a(i,j,k,1) = 0.0;//sin( 2*(geom.ProbLo(X) + i*geom.CellSize(X) +1)*3.1415962);
         a(i,j,k,2) = 0.5; 
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
    // Simulation parameters,  these should be read from a file quite soon
    
    const  int n_cell = 32;
    int max_grid_size = 16;
    int nsteps = 3000;
    double dt = 1.0/250;
    // Do a quite even load balancing
    amrex::DistributionMapping::strategy(amrex::DistributionMapping::KNAPSACK);

    // Periodic
    amrex::Vector<int> is_periodic(AMREX_SPACEDIM,1);     
    // Nodal indexing
    amrex::IndexType typ({AMREX_D_DECL(0,0,0)});


    amrex::IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
    amrex::IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, 32-1));
    amrex::Box domain(dom_lo, dom_hi,typ);
    amrex::BoxArray ba(domain);

    // Initialize the boxarray "ba" from the single box "bx"
    // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
    ba.maxSize(max_grid_size);

    // This defines the physical box, [-1,1] in each direction.
    amrex::RealBox real_box({AMREX_D_DECL(-1,-1,-1)},
                     {AMREX_D_DECL(1,1,1)});

    // This defines a Geometry object
    amrex::Geometry geom(domain,&real_box,amrex::CoordSys::cartesian,is_periodic.data());
    // How Boxes are distrubuted among MPI processes
    amrex::DistributionMapping dm(ba);

    int Nghost = 2;
    
    // Ncomp = number of components for each array
    int Ncomp  = 3;
  
    amrex::MultiFab E(ba,dm,Ncomp,Nghost);
    amrex::MultiFab B(ba,dm,Ncomp,Nghost);
    CParticleContainer P(geom,dm,ba,3);

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

    
for(amrex::MFIter mfi= P.MakeMFIter(0) ;mfi.isValid();++mfi){
    
    // Each grid,tile has a their own local particle container
    auto& particles = P.GetParticles(0)[std::make_pair(mfi.index(),
                                        mfi.LocalTileIndex())];
    auto box=mfi.validbox();
    const auto lo = amrex::lbound(box);
    //const auto hi = amrex::ubound(box);
    if(mfi.index()==0){
    add_single_particle(particles,{0.0,-0.2,0},{0.1,0,0},1,-1);
//    add_single_particle(particles,{-0.9,-0.9,-0.9},{0.1,1,1},10,-1);
//    add_single_particle(particles,{-0.95,-0.95,-0.95},{0.1,1,1},10,-1);
    }
    //  add_single_particle(particles,{0.001,0,0},{0.1,0,0},10,-1);
}

    P.Redistribute();
    P.fillNeighbors();
    P.updateNeighbors();

    int id = amrex::ParallelDescriptor::MyProc();
    print_Particle_info(geom,P);


for(int step=0; step<nsteps;step++){


    if(id == 0){
    
    std::cout <<"Step:" <<step << std::endl;
    } 

    auto E_tot = get_total_energy(geom,P,E,B); 
   amrex::Print() <<"ENERGY: "<<E_tot.first <<" "<< E_tot.second << std::endl;
/*    if(step % 5 ==0){

    int n=step;
    amrex::Real time=step*dt;
    const std::string& pltfile_E = amrex::Concatenate("plt_E",n,0);
    WriteSingleLevelPlotfile(pltfile_E, E, {"E_x","E_y","E_z"}, geom, time, n);
    const std::string& pltfile_B = amrex::Concatenate("plt_B",n,0);
    WriteSingleLevelPlotfile(pltfile_B, B, {"B_x","B_y","B_z"}, geom, time, n);
  
    }
*/
    G_Theta_E(geom,P,E,B,dt/2);
    G_Theta<X>(geom,P,E,B,dt/2);
    G_Theta<Y>(geom,P,E,B,dt/2);
    G_Theta<Z>(geom,P,E,B,dt/2);
//    G_Theta_B(geom,P,E,B,dt); 
    G_Theta<Z>(geom,P,E,B,dt/2);
    G_Theta<Y>(geom,P,E,B,dt/2);
    G_Theta<X>(geom,P,E,B,dt/2);
    G_Theta_E(geom,P,E,B,dt/2);
    print_Particle_info(geom,P);



}

//    print_Particle_info(geom,P);
//    auto E_tot = get_total_energy(geom,P,E,B); 
//    std::cout <<"ENERGY: "<<E_tot.first <<" "<< E_tot.second << std::endl;


}

