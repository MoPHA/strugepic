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
#include <AMReX_BC_TYPES.H>
#include <AMReX_BCRec.H>
#include <AMReX_BCUtil.H>


// std c++
#include <iostream>
#include <malloc.h>
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
        if(i > 199 &&i < 251){
         a(i,j,k,0) = 0;
         a(i,j,k,1) = 0.25*sin( 2*((1.0*i)/50)*M_PI); 
         a(i,j,k,2) = 0; 
        }
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
        if(i > 199 &&i < 251){
         a(i,j,k,0) = 0;
         a(i,j,k,1) = 0;
         a(i,j,k,2) = -0.25*sin( 2*((1.0*i)/50)*M_PI); 
        }
       }
     }
   }
}

void DebugPrint(const amrex::Geometry geom, amrex::MultiFab &A){
   
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

                std::cout << a(i,j,k,Y) << std::endl;
            }
         });

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
    
    const  int n_cell = 300;
    int max_grid_size = 300;
    int nsteps = 1000;
    double dt = 0.5;
    double q=-4.80320467059932e-11;
    double m=1.5453871347313696e-07;
    double v=0.020013845711889123;
    
   
    double omega = 2*M_PI/120;
    double source_pos=4;
    
    auto Es = E_source(source_pos,omega,dt);

    // Do a quite even load balancing
    amrex::DistributionMapping::strategy(amrex::DistributionMapping::KNAPSACK);

    // Periodic
    amrex::Vector<int> is_periodic({0,1,1});     


    
    // cell centered indexing
    // We acctually calculate as if the fields are centered at the nodes
    amrex::IndexType typ({AMREX_D_DECL(0,0,0)});


    amrex::IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
    amrex::IntVect dom_hi(AMREX_D_DECL(n_cell-1, 5-1, 5-1));
    amrex::Box domain(dom_lo, dom_hi,typ);
    amrex::BoxArray ba(domain);

    // Initialize the boxarray "ba" from the single box "bx"
    // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
    ba.maxSize(max_grid_size);

    // This defines the physical box, [-1,1] in each direction.
    amrex::RealBox real_box({AMREX_D_DECL(0,0,0)},
                     {AMREX_D_DECL(n_cell,5,5)});

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
//    add_single_particle(particles,{0.18194016191876186,0.01,0.01},{v,0,0},m,q);
//    add_single_particle(particles,{-0.9,-0.9,-0.9},{0.1,1,1},10,-1);
//    add_single_particle(particles,{-0.95,-0.95,-0.95},{0.1,1,1},10,-1);
    }
    //  add_single_particle(particles,{0.001,0,0},{0.1,0,0},10,-1);
}

//    double m=1;
//    double q=-1;
//    add_particle_one_per_cell(geom,P,m,q);


    P.Redistribute();
    P.fillNeighbors();
    P.updateNeighbors();

   int id = amrex::ParallelDescriptor::MyProc();
for(int step=0; step<nsteps;step++){
   // wave_source(geom,E,B,step*dt);

    if(id == 0){
    
    std::cout <<"Step:" <<step << std::endl;
    } 

    auto E_tot = get_total_energy(geom,P,E,B); 
    amrex::Print() <<"ENERGY: "<<E_tot.first <<" "<< E_tot.second << std::endl;
    if(step % 10 ==0){

    int n=step;
    amrex::Real time=step*dt;
    const std::string& pltfile_E = amrex::Concatenate("plt_E",n,0);
    WriteSingleLevelPlotfile(pltfile_E, E, {"E_x","E_y","E_z"}, geom, time, n);
    const std::string& pltfile_B = amrex::Concatenate("plt_B",n,0);
    WriteSingleLevelPlotfile(pltfile_B, B, {"B_x","B_y","B_z"}, geom, time, n);
    P.WriteBinaryParticleData(amrex::Concatenate("plt_P",step,0),"Particle0",{1,1,1,1,1},{},{"Mass","Charge","VX","VY","VZ"},{});
    }

    G_Theta_E(geom,P,E,B,dt/2);
    G_Theta<X>(geom,P,E,B,dt/2);
    G_Theta<Y>(geom,P,E,B,dt/2);
    G_Theta<Z>(geom,P,E,B,dt/2);
    Es(geom,E,step*dt);
    G_Theta_B(geom,P,E,B,dt);
    G_Theta<Z>(geom,P,E,B,dt/2);
    G_Theta<Y>(geom,P,E,B,dt/2);
    G_Theta<X>(geom,P,E,B,dt/2);
    G_Theta_E(geom,P,E,B,dt/2);
    print_Particle_info(geom,P);



}


}

