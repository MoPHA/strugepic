#include "AMReX_Array.H"
#include "AMReX_Box.H"
#include "AMReX_BoxArray.H"
#include "AMReX_DistributionMapping.H"
#include "AMReX_Geometry.H"
#include "AMReX_MultiFab.H"
#include "AMReX_ParGDB.H"
#include "AMReX_ParallelDescriptor.H"
#include <AMReX_PlotFileUtil.H>
#include "AMReX_ParallelReduce.H"
#include "AMReX_REAL.H"
#include <iostream>
#include <math.h>
#include <utility>
#include <vector>
#include "AMReX_VisMF.H"
#include "particle_defs.hpp"
#include "amrex_util.hpp"
#include <fstream>
#include "w_defs.hpp"
#include "propagators.hpp"

// Init external dirichle boundary condition

// Either periodic or dirichle

void set_uniform_field(amrex::MultiFab &A, std::array<double,3> vals){

    for (amrex::MFIter mfi(A); mfi.isValid(); ++mfi){
        const amrex::Box& box = mfi.validbox();
        amrex::Array4<amrex::Real> const& a = A.array(mfi); 

    amrex::ParallelFor(box,  [=] AMREX_GPU_DEVICE (int i,int j,int k ){

                a(i,j,k,X)= vals[X];
                a(i,j,k,Y)= vals[Y];
                a(i,j,k,Z)= vals[Z];

            });

    }
}

void SimulationIO::dump_pdens(std::string filename){

    get_particle_number_density<WRANGE>(geom,P,Pdens); 
    for (amrex::MFIter mfi(Pdens); mfi.isValid(); ++mfi){
        const amrex::Box& box = mfi.validbox();
        amrex::Array4<amrex::Real> const& pd = Pdens.array(mfi); 

    amrex::ParallelFor(box,  [=] AMREX_GPU_DEVICE (int i,int j,int k ){

            amrex::AllPrintToFile(filename) <<i << " "<< j << " " << k << " "  <<  pd(i,j,k,0) << "\n";

            });

    }

}

inline void dump_field(amrex::MultiFab & A , std::string filename){
    for (amrex::MFIter mfi(A); mfi.isValid(); ++mfi){
        const amrex::Box& box = mfi.validbox();
        amrex::Array4<amrex::Real> const& a = A.array(mfi); 
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i,int j,int k ){
            amrex::AllPrintToFile(filename) <<i << " "<< j << " " << k << " "  <<  a(i,j,k,X) << " " <<a(i,j,k,Y) << " " << a(i,j,k,Z) << "\n";
            });

    }

}

void SimulationIO::dump_E_field(std::string filename){
        dump_field(this->E,filename);   
}
void SimulationIO::dump_B_field(std::string filename){
        dump_field(this->B,filename);   
}



void SimulationIO::write(int step,bool checkpoint,bool particles){
    if(checkpoint){
    amrex::VisMF::Write(E,amrex::Concatenate(data_folder_name+std::string("/E_CP"),step,0));
    amrex::VisMF::Write(B,amrex::Concatenate(data_folder_name+std::string("/B_CP"),step,0));
    P.Checkpoint(amrex::Concatenate(data_folder_name+std::string("/P_CP"),step,0),"Particle0");
    }
    else{
    get_particle_number_density<WRANGE>(geom,P,Pdens);
    int n=step;
    amrex::Real time=step*dt;
    const std::string& pltfile_E = amrex::Concatenate(data_folder_name+std::string("/plt_E"),n,0);
    WriteSingleLevelPlotfile(pltfile_E, E, {"E_x","E_y","E_z"}, geom, time, n);
    const std::string& pltfile_B = amrex::Concatenate(data_folder_name+std::string("/plt_B"),n,0);
    WriteSingleLevelPlotfile(pltfile_B, B, {"B_x","B_y","B_z"}, geom, time, n);
    const std::string& pltfile_Pdens = amrex::Concatenate(data_folder_name+std::string("/plt_Pdens"),n,0);
    WriteSingleLevelPlotfile(pltfile_Pdens, Pdens, {"n"}, geom, time, n);
    }
    if(particles){
        P.WriteBinaryParticleData(amrex::Concatenate(data_folder_name+std::string("/plt_P"),step,0),"Particle0",{1,1,1,1,1},{},{"Mass","Charge","VX","VY","VZ"},{});
    }
}
void SimulationIO::read(int step){
    amrex::VisMF::Read(E,amrex::Concatenate(data_folder_name+std::string("/E_CP"),step,0));
    amrex::VisMF::Read(B,amrex::Concatenate(data_folder_name+std::string("/B_CP"),step,0));
    P.Restart(amrex::Concatenate(data_folder_name+std::string("/P_CP"),step,0),"Particle0");
}



SimulationIO::SimulationIO(amrex::Geometry geom,amrex::MultiFab & E,amrex::MultiFab & B,CParticleContainer &P,double dt,std::string data_folder_name):
    geom(geom),E(E),B(B),P(P),dt(dt),data_folder_name(data_folder_name){
        int Nghost=E.nGrow();
        this->Pdens=amrex::MultiFab(P.ParticleBoxArray(0),P.ParticleDistributionMap(0),1,Nghost);
    }

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

void add_single_particle(CParticleContainer&P,amrex::RealArray pos , amrex::RealArray vel, double m,double q){
for(amrex::MFIter mfi= P.MakeMFIter(0) ;mfi.isValid();++mfi){
    
    // Each grid,tile has a their own local particle container
    auto& particles = P.GetParticles(0)[std::make_pair(mfi.index(),mfi.tileIndex())];
        if(mfi.index()==0){
            add_single_particle(particles,pos,vel,m,q);
        }
    }
    P.Redistribute();

}


// The dist function should return [0,1]
// m and q are the normalized mass and charge for the "real" particle
// These are then scaled according to the number of computational particles 

double bernstein_density(const amrex::Geometry geom, int i , int j, int k){
    int nr=380;
    if(i < nr+320){
        return exp( -(i-(nr+320))*(i-(nr+320))/( 2*(nr/3.5)*(nr/3.5) ) ); 
    }
    else if(i >= 1300){
          return exp(-(i-1300)*(i-1300)/26122.0);
    }
    else{
    return 1;
    }

}
double uniform_density(const amrex::Geometry geom,int i ,int j ,int k){
    return 1;
}

double simple_line_density(const amrex::Geometry geom ,int i , int j , int k){
      return (1.0*i/20);  
}

double gaussian_dist(double pos,double center,double std_dev){
        // Scale and shift position to center
        double sc_pos=(pos-center)/(std_dev);
        return exp(-sc_pos*sc_pos);
}

double d_gaussian_dist(double pos,double center,double std_dev){
            return -2*(pos-center)/(std_dev*std_dev)*gaussian_dist(pos,center,std_dev);
}

void set_field_gradient_gaussian_x(amrex::MultiFab &A,double A_max,double center,double std_dev){

    for (amrex::MFIter mfi(A); mfi.isValid(); ++mfi){
        const amrex::Box& box = mfi.validbox();
        amrex::Array4<amrex::Real> const& a = A.array(mfi); 

    amrex::ParallelFor(box,  [=] AMREX_GPU_DEVICE (int i,int j,int k ){
                a(i,j,k,X)=A_max*d_gaussian_dist(i,center,std_dev);                 
            });

    }

}



// This could be distributed??
void distribute_processes_pdens(amrex::DistributionMapping dm,const amrex::Geometry geom,amrex::BoxArray &Ba,double (*dist_func)(const amrex::Geometry,int,int,int),std::string strat){
        int idx=0;
        std::vector<long int> I(Ba.boxList().size());
    for(auto b:Ba.boxList()){
        auto lo=amrex::lbound(b);
        auto hi=amrex::ubound(b);
        for(int k= lo.z ;k<hi.z; k++){
        for(int j= lo.y ;j<hi.y; j++){
        for(int i= lo.x ;i<hi.x ; i++){
                I[idx]+=dist_func(geom,i,j,k)*10000;
        }
        }
        }

    idx+=1;
    }
    if(strat=="SFC"){
    dm.SFCProcessorMap(Ba,I,amrex::ParallelDescriptor::NProcs());
    }
    else if(strat=="KnapSack"){
    dm.KnapSackProcessorMap(I,amrex::ParallelDescriptor::NProcs());
    }
    else{
    amrex::Print() << "Unknown strategy!!, supported are SFC and KnapSack";
    exit(1);
    }

}


void add_particle_density(const amrex::Geometry geom , CParticleContainer&P, double (*dist_func)(const amrex::Geometry,int,int,int),int ppc_max ,double m, double q, double v){

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0,1);
    std::normal_distribution<double> vel(0,v); 
    double q_c = q/ppc_max;
    double m_c = m/ppc_max;

// For simplicity, particles are initialized with Maxwell–Boltzmann distribution


for(amrex::MFIter mfi= P.MakeMFIter(0) ;mfi.isValid();++mfi){
    
    // Each grid,tile has a their own local particle container
    auto& particles = P.GetParticles(0)[std::make_pair(mfi.index(),
                                        mfi.LocalTileIndex())];
    auto box=mfi.validbox();

    
   auto domain=geom.Domain();
   auto lod=amrex::lbound(domain);
   auto hid=amrex::ubound(domain);


   const auto lo = amrex::lbound(box);
   const auto hi = amrex::ubound(box);
   for     (int k = lo.z; k <= hi.z; ++k) {
     for   (int j = lo.y; j <= hi.y; ++j) {
       for (int i = lo.x; i <= hi.x; ++i) { 
            
           double x = geom.ProbLo(X) + i*geom.CellSize(X);
           double y = geom.ProbLo(Y) + j*geom.CellSize(Y);
           double z = geom.ProbLo(Z) + k*geom.CellSize(Z);
           int num_particles = dist_func(geom,i,j,k)*ppc_max;  
           if((i <= PART_BOUND+lod.x  || i >= hid.x-PART_BOUND) && !geom.isPeriodic(X)){
            continue;
           }
            for(int p =0; p < num_particles ; p++){ 
                add_single_particle(particles,{x+dist(mt)*geom.CellSize(X),y+dist(mt)*geom.CellSize(Y),z+dist(mt)*geom.CellSize(Z)},{vel(mt),vel(mt),vel(mt)},m_c,q_c);
            }
       }
     }
   }
 }
P.Redistribute();

}


void add_particle_n_per_cell(amrex::Geometry geom, CParticleContainer&P,double m,double q,double v,int n){
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0,1);
    std::normal_distribution<double> vel(0,v); 
for(amrex::MFIter mfi= P.MakeMFIter(0) ;mfi.isValid();++mfi){
    
    // Each grid,tile has a their own local particle container
    auto& particles = P.GetParticles(0)[std::make_pair(mfi.index(),
                                        mfi.LocalTileIndex())];
    auto box=mfi.validbox();

   auto domain=geom.Domain();
   auto lod=amrex::lbound(domain);
   auto hid=amrex::ubound(domain);

   const auto lo = amrex::lbound(box);
   const auto hi = amrex::ubound(box);
   for     (int k = lo.z; k <= hi.z; ++k) {
     for   (int j = lo.y; j <= hi.y; ++j) {
       for (int i = lo.x; i <= hi.x; ++i) { 
            
           double x = geom.ProbLo(X) + i*geom.CellSize(X);
           double y = geom.ProbLo(Y) + j*geom.CellSize(Y);
           double z = geom.ProbLo(Z) + k*geom.CellSize(Z);
            for(int p =0; p < n ; p++){ 
           if((i <= PART_BOUND+lod.x  || i >= hid.x-PART_BOUND) && !geom.isPeriodic(X)){
            continue;
           }
                add_single_particle(particles,{x+dist(mt)*geom.CellSize(X),y+dist(mt)*geom.CellSize(Y),z+dist(mt)*geom.CellSize(Z)},{vel(mt),vel(mt),vel(mt)},m/n,q/n);
            }
       }
     }
   }
 }
P.Redistribute();

}


void print_Particle_info(const amrex::Geometry geom,CParticleContainer&P ){

    for (CParIter pti(P, 0); pti.isValid(); ++pti) {
        auto&  particles = pti.GetArrayOfStructs();
        for(auto p : particles ){
            std::cout << "(" << p.cpu()<<"," << p.id()<<"," << amrex::ParallelDescriptor::MyProc() <<")" << std::endl;
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
