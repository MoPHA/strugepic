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

void main_main();

void init_E (amrex::Box const& bx, amrex::Array4<amrex::Real> const& a)
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
void init_B (amrex::Box const& bx, amrex::Array4<amrex::Real> const& a)
{
   const auto lo = amrex::lbound(bx);
   const auto hi = amrex::ubound(bx);
   for     (int k = lo.z; k <= hi.z; ++k) {
     for   (int j = lo.y; j <= hi.y; ++j) {
       for (int i = lo.x; i <= hi.x; ++i) { 
         a(i,j,k,0) = 0;
         a(i,j,k,1) = 0;
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
    
    int n_cell = 64;
    int max_grid_size=32;
    int nsteps=100;
    // Do a quite even load balancing
    amrex::DistributionMapping::strategy(amrex::DistributionMapping::KNAPSACK);

    // Periodic
    amrex::Vector<int> is_periodic(AMREX_SPACEDIM,1);     
    // Nodal indexing
    amrex::IndexType typ({AMREX_D_DECL(1,1,1)});


    amrex::IntVect dom_lo(AMREX_D_DECL(       0,        0,        0));
    amrex::IntVect dom_hi(AMREX_D_DECL(n_cell-1, n_cell-1, n_cell-1));
    amrex::Box domain(dom_lo, dom_hi,typ);
    amrex::BoxArray ba(domain);

    // Initialize the boxarray "ba" from the single box "bx"
    // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
    ba.maxSize(max_grid_size);

    // This defines the physical box, [-1,1] in each direction.
    amrex::RealBox real_box({AMREX_D_DECL(0,0,0)},
                     {AMREX_D_DECL( n_cell,n_cell,n_cell)});

    // This defines a Geometry object
    amrex::Geometry geom(domain,&real_box,amrex::CoordSys::cartesian,is_periodic.data());
    
    // How Boxes are distrubuted among MPI processes
    amrex::DistributionMapping dm(ba);

    int Nghost = 2;
    
    // Ncomp = number of components for each array
    int Ncomp  = 3;
  
    amrex::MultiFab E(ba,dm,Ncomp,Nghost);
    amrex::MultiFab B(ba,dm,Ncomp,Nghost);
    amrex::NeighborParticleContainer<5,0> P(geom,dm,ba,3);

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
    init_E(box, a);
    init_B(box, b);
}

    // Init single particle
    
for(amrex::MFIter mfi= P.MakeMFIter(0) ;mfi.isValid();++mfi){
    
    // Each grid,tile has a their own local particle container
    auto& particles = P.GetParticles(0)[std::make_pair(mfi.index(),
                                        mfi.LocalTileIndex())];
    //if(mfi.index() !=0){
    //    continue;
   // }
   
    auto box=mfi.validbox();
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    amrex::Particle<5> p;
    p.id()   = amrex::Particle<5>::NextID();
    p.cpu()  = amrex::ParallelDescriptor::MyProc();
    p.pos(0) = lo.x+1.1;
    p.pos(1) = lo.y+1.1;
    p.pos(2) = lo.z+1.1;
    p.rdata(0)=1;
    p.rdata(1)=-1;
    p.rdata(2)=0.5;
    p.rdata(3)=0;
    p.rdata(4)=0;
    particles.push_back(p);
}

    P.Redistribute();
    P.fillNeighbors();
    P.updateNeighbors();

    E.FillBoundary(geom.periodicity());
    B.FillBoundary(geom.periodicity());

using MyParIter = amrex::ParConstIter<5,0,0,0>;
for (MyParIter pti(P, 0); pti.isValid(); ++pti) {
    const auto& particles = pti.GetArrayOfStructs();
    amrex::FArrayBox& efab = E[pti];
    const amrex::Box& box = pti.validbox();;
    amrex::Array4<amrex::Real> const& a = efab.array();
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    amrex::Print() << pti.validbox() << std::endl;
    for (const auto& p : particles) {
    amrex::Print() << p.pos(0) << " " << p.pos(1) << " " << p.pos(2) << std::endl ;
    }
    const auto& n_particles = P.GetNeighbors(0,pti.index(),pti.LocalTileIndex());
    for(const auto& p: n_particles){
    }
}

    





    int n=0;
    amrex::Real time=0.0;
    const std::string& pltfile_E = amrex::Concatenate("plt_E",n,2);
    WriteSingleLevelPlotfile(pltfile_E, E, {"E_x","E_y","E_z"}, geom, time, n);
    const std::string& pltfile_B = amrex::Concatenate("plt_B",n,2);
    WriteSingleLevelPlotfile(pltfile_B, B, {"B_x","B_y","B_z"}, geom, time, n);
    const std::string& pltfile_P = amrex::Concatenate("plt_P",n,2);
    
    const amrex::Vector<std::string> var_names={"Mass","Charge","Vx","Vy","Vz"};
    const std::string p_name="test_particle";
    P.WritePlotFile(pltfile_P,p_name,var_names);



}

