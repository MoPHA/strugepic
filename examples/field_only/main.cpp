// Amrex
#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_RealBox.H>
#include <AMReX_CoordSys.H>
#include <AMReX_MultiFab.H>
#include <AMReX_DistributionMapping.H>
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
#include "strugepic_util.hpp"
#include "strugepic_propagators.hpp"
#include "strugepic_defs.hpp"
#include "cmath"

void set_sin_field(amrex::MultiFab &A){

    for (amrex::MFIter mfi(A); mfi.isValid(); ++mfi){
        const amrex::Box& box = mfi.validbox();
        amrex::Array4<amrex::Real> const& a = A.array(mfi); 

    amrex::ParallelFor(box,  [=] AMREX_GPU_DEVICE (int i,int j,int k ){

                a(i,j,k,Y)= sin(2*i*2*3.14151962/300);

            });

    }
}

void main_main();

int main(int argc, char* argv[])
{

    amrex::Initialize(argc,argv);
    main_main();
    amrex::Finalize();

}
void main_main()
{
    amrex::ParmParse pp;

    std::array<int,3> n_cell;
    std::array<int,3> max_grid_size;
    int x_periodic;
    int Nghost = 3; 

    int nsteps;
    int start_step;
    double dt;
    int output_interval;
    int checkpoint_interval;
    
    double Es;
    double omega;
    int sp;
    int Ncomp  = 3;

    std::string data_folder_name;


    pp.get("output_interval",output_interval);
    pp.get("checkpoint_interval",checkpoint_interval);
    pp.get("n_cell",n_cell);
    pp.get("max_grid_size",max_grid_size);
    pp.get("x_periodic",x_periodic);
    pp.get("nsteps",nsteps);
    pp.get("start_step",start_step);
    pp.get("sp",sp);
    pp.get("dt",dt);
    pp.get("Es",Es);
    pp.get("omega",omega);
    pp.get("data_folder_name",data_folder_name);







    


///////////7



    // Periodic
    amrex::Vector<int> is_periodic({x_periodic,1,1});     
    // Cell indexing
    amrex::IndexType typ({AMREX_D_DECL(0,0,0)});


    amrex::IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
    amrex::IntVect dom_hi(AMREX_D_DECL(n_cell[X]-1, n_cell[Y]-1, n_cell[Z]-1));
    amrex::Box domain(dom_lo, dom_hi,typ);
    amrex::BoxArray ba(domain);

    // Initialize the boxarray "ba" from the single box "bx"
    // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
    ba.maxSize({max_grid_size[X],max_grid_size[Y],max_grid_size[Z]});
    // This defines the physical box, [-1,1] in each direction.
    amrex::RealBox real_box({AMREX_D_DECL(0,0,0)},
                     {AMREX_D_DECL((double)n_cell[X] , (double)n_cell[Y],(double)n_cell[Z])});

    // This defines a Geometry object
    amrex::Geometry geom(domain,&real_box,amrex::CoordSys::cartesian,is_periodic.data());
    // How Boxes are distrubuted among MPI processes
    amrex::DistributionMapping dm(ba);
    CParticleContainer P(geom,dm,ba);
    //distribute_processes_pdens(dm,geom,ba,bernstein_density,"SFC");    


    amrex::MultiFab E(ba,dm,Ncomp,Nghost);
    amrex::MultiFab B(ba,dm,Ncomp,Nghost);
    auto SimIO=SimulationIO(geom,E,B,P,dt,data_folder_name);
    auto Source=E_source(geom,E,sp,Y,Es,omega,dt);


////////////7

    if(start_step !=0){
    SimIO.read(start_step);
    }

    E.FillBoundary(geom.periodicity());
    B.FillBoundary(geom.periodicity());


for(int step=start_step; step<nsteps;step++){ 
    print_Particle_info(geom,P);
    amrex::Print() <<"Step:" <<step << std::endl;
    auto E_tot = get_total_energy(geom,P,E,B); 
    amrex::Print() <<"ENERGY: "<<E_tot.first <<" "<< E_tot.second << std::endl;
    if(step % output_interval ==0 && output_interval != -1){
        SimIO.write<WRANGE>(step);
    }
    if(step % checkpoint_interval ==0 && checkpoint_interval  !=-1){
        SimIO.write<WRANGE>(step,true,false);
    }

    G_Theta_E<WRANGE>(geom,P,E,B,dt/2);
    Source(dt*step);
    G_Theta_B(geom,P,E,B,dt);
    G_Theta_E<WRANGE>(geom,P,E,B,dt/2);
}



}

