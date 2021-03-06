#ifndef UTIL
#define UTIL
#include "AMReX_Array.H"
#include "AMReX_Geometry.H"
#include "AMReX_MultiFab.H"
#include <AMReX_PlotFileUtil.H>
#include "strugepic_defs.hpp"
#include "strugepic_w.hpp"
#include <tuple>
#include <AMReX_GpuQualifiers.H>
#include <math.h>
// What cell index is a given point in?
// This is equivalent with the index for the "Lower left" corner
amrex::IntArray get_point_cell(const amrex::Geometry geom,const amrex::RealArray pos);
void add_single_particle( CParticleTile&particlet ,amrex::RealArray pos , amrex::RealArray vel, double m,double q);
void FillDirichletBoundary(const amrex::Geometry geom, amrex::MultiFab &A );
double bernstein_density(const amrex::Geometry geom, int i , int j, int k);
double simple_line_density(const amrex::Geometry geom ,int i , int j , int k);
void add_particle_density(const amrex::Geometry geom , CParticleContainer&P, double (*dist_func)(const amrex::Geometry,int,int,int),int ppc_max ,double m, double q ,double v);
void distribute_processes_pdens(amrex::DistributionMapping dm,const amrex::Geometry geom,amrex::BoxArray &Ba,double (*dist_func)(const amrex::Geometry,int,int,int),std::string strat);
void set_uniform_field(amrex::MultiFab &A, std::array<double,3> vals);
void add_single_particle(CParticleContainer&P,amrex::RealArray pos , amrex::RealArray vel, double m,double q);
void print_boxes(amrex::BoxArray ba);
double uniform_density(const amrex::Geometry geom,int i ,int j ,int k);
double gaussian_dist(double pos,double center,double std_dev);
void set_field_gradient_gaussian_x(amrex::MultiFab &A,double A_max,double center,double std_dev);

AMREX_GPU_HOST_DEVICE int construct_segments(amrex::Real x_start ,amrex::Real x_end, amrex::Real *seg_points,int *seg_idx);

template<int W_range>
void get_particle_number_density(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab &P_dens){
    P_dens.setVal(0);
    P_dens.setBndry(0);
    for (CParIter pti(P,0); pti.isValid(); ++pti){

        // Each grid,tile has a their own local particle container
        CParticle * AMREX_RESTRICT particles=&(pti.GetArrayOfStructs()[0]);
        const int np= pti.numParticles();
        amrex::Array4<amrex::Real> const& P_dens_loc = P_dens[pti].array();


    const auto _Ics = geom.InvCellSize();
    // Remote particle grid lower corner
    amrex::GpuArray<amrex::Real,3> lb;
    amrex::GpuArray<amrex::Real,3> Ics;
    lb[X] =geom.ProbLo(X);
    lb[Y] =geom.ProbLo(Y);
    lb[Z] =geom.ProbLo(Z);
    Ics[X]=_Ics[X];
    Ics[Y]=_Ics[Y];
    Ics[Z]=_Ics[Z];

    amrex::ParallelFor(np,[=] AMREX_GPU_DEVICE (long p_i) {
        int coord[3];
        coord[X]=floor((particles[p_i].pos(X) -lb[X])*Ics[X]);
        coord[Y]=floor((particles[p_i].pos(Y) -lb[Y])*Ics[Y]);
        coord[Z]=floor((particles[p_i].pos(Z) -lb[Z])*Ics[Z]);



        auto px=(particles[p_i].pos(X)-lb[X])*Ics[X];
        auto py=(particles[p_i].pos(Y)-lb[Y])*Ics[Y];
        auto pz=(particles[p_i].pos(Z)-lb[Z])*Ics[Z];

        constexpr int W1_li=-W_range+1;
        constexpr int W1_hi= W_range;
        constexpr int Wp_li= W1_li;
        constexpr int Wp_hi= W1_hi-1;

        for(int i=Wp_li;i<=Wp_hi;i++){
            auto nx=coord[X] +i;
            for(int j=Wp_li;j<=Wp_hi;j++){
                auto ny=coord[Y] +j;
                for(int k=Wp_li;k<=Wp_hi;k++){
                    auto nz=coord[Z] +k;
                       amrex::Gpu::Atomic::Add(&P_dens_loc(nx,ny,nz,0),Wp(px-nx)*Wp(py-ny)*Wp(pz-nz));
            }
        }
        }

    });

    }
    P_dens.SumBoundary(geom.periodicity());
}

template<int comp>
void shift_periodic(const amrex::Geometry geom ,CParticle &particle){

    double pos=particle.pos(comp);
    double lb =geom.ProbLo(comp);
    double ub =geom.ProbHi(comp);

    if(pos > ub){
        particle.pos(comp) = lb+(pos-ub);
    }
    else if(pos < lb){

        particle.pos(comp) = ub-(lb-pos);
    }


}

class SimulationIO
{
    public:
    SimulationIO(amrex::Geometry geom,amrex::MultiFab & E,amrex::MultiFab & B,CParticleContainer &P,double dt,std::string data_folder_name);
    template <int W_range>
    void write(int step,bool checkpoint=false,bool particles=false);
    void read(int step);
    void dump_E_field(std::string filename);
    void dump_B_field(std::string filename);
    private:
        amrex::Geometry geom;
        amrex::MultiFab & E;
        amrex::MultiFab & B;
        amrex::MultiFab  Pdens;
        CParticleContainer &P;
        double dt;
        std::string data_folder_name;


};

template <int W_range>
void SimulationIO::write(int step,bool checkpoint,bool particles){
    if(checkpoint){
    amrex::VisMF::Write(E,amrex::Concatenate(data_folder_name+std::string("/E_CP"),step,0));
    amrex::VisMF::Write(B,amrex::Concatenate(data_folder_name+std::string("/B_CP"),step,0));
    P.Checkpoint(amrex::Concatenate(data_folder_name+std::string("/P_CP"),step,0),"Particle0");
    }
    else{
    get_particle_number_density<W_range>(geom,P,Pdens);
    int n=step;
    amrex::Real time=step*dt;
    const std::string& pltfile_E = amrex::Concatenate(data_folder_name+std::string("/plt_E"),n,0);
    amrex::WriteSingleLevelPlotfile(pltfile_E, E, {"E_x","E_y","E_z"}, geom, time, n);
    const std::string& pltfile_B = amrex::Concatenate(data_folder_name+std::string("/plt_B"),n,0);
    amrex::WriteSingleLevelPlotfile(pltfile_B, B, {"B_x","B_y","B_z"}, geom, time, n);
    const std::string& pltfile_Pdens = amrex::Concatenate(data_folder_name+std::string("/plt_Pdens"),n,0);
    amrex::WriteSingleLevelPlotfile(pltfile_Pdens, Pdens, {"n"}, geom, time, n);
    }
    if(particles){
        amrex::Print() << "Writing binary particle data not implemented due to change in amrex api" << std::endl;
    }
}




template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

void print_Particle_info(const amrex::Geometry geom,CParticleContainer&P );



template<int comp>
int get_point_line(const amrex::Geometry geom,const amrex::Real pos){
    const auto icellsize=geom.InvCellSize(comp);
    auto problo=geom.ProbLo(comp);
    return floor((pos -problo)*icellsize);
}






template<int comp,int W_range>
AMREX_GPU_HOST_DEVICE inline bool segment_reflect(int isPeriodic ,int smallEnd, int bigEnd,amrex::Real * seg_points,int * seg_idx){
if(!isPeriodic && ((seg_idx[1] == smallEnd+W_range ) || (seg_idx[1] == bigEnd-W_range) )){
        seg_idx[1]=seg_idx[0];
        seg_points[2]=2*seg_points[1]-seg_points[2];
        return true;
    }
        return false;
}

template<int comp>
AMREX_GPU_HOST_DEVICE inline void particle_reflect(CParticle *p,amrex::Real* seg_points){
    p->pos(comp)=seg_points[2];
    p->rdata(comp+2)*=-1;
}




std::pair<amrex::Real,amrex::Real> get_total_energy(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab &E, amrex::MultiFab &B );
void add_particle_n_per_cell(const amrex::Geometry geom, CParticleContainer&P,double m,double q,double v,int n);
template <int comp>
void set_two_stream_drift(CParticleContainer &P,double v){
    for(CParIter pti(P,0);pti.isValid() ; ++pti){
        auto & particles = pti.GetArrayOfStructs();
        for (auto &p : particles){
             p.rdata(comp +2)+= ((p.id() % 2)*2 -1)*v;
        }
    }

}

template <int comp>
void set_temprature_gradient_gaussian(CParticleContainer &P,double v_min,double v_max,double center,double std_dev,double B_init){
    std::random_device rd;
    std::mt19937 mt(rd());

    double v_min_var=v_min*v_min;
    double v_max_var=v_max*v_max;
    double v_diff_var=v_max_var-v_min_var;
    std::uniform_real_distribution<double> angle(0,2*M_PI);
    for(CParIter pti(P,0);pti.isValid() ; ++pti){
        auto & particles = pti.GetArrayOfStructs();
        for (auto &p : particles){
            double v = sqrt(v_min_var+gaussian_dist(p.pos(comp),center,std_dev)*v_diff_var);
            double ang = angle(mt);
            std::normal_distribution<double> vel(0,v);
            double vx = vel(mt);
            double vy = vel(mt);
            double v_mag=sqrt(vx*vx+vy*vy);
            double rad = p.rdata(M)*v_mag/(p.rdata(Q)*B_init);

            p.rdata(VX)= -v_mag*sin(ang);
            p.rdata(VY)= v_mag*cos(ang);
            p.pos(X)+=rad*cos(ang);
            p.pos(Y)+=rad*sin(ang);

            p.rdata(VZ)= vel(mt);
        }
    }

}

template  <int p_comp,int v_comp>
void save_velocity_distribution( CParticleContainer &P, double p_min,double p_max, double v_min,double v_max,int p_bins,int v_bins,std::string filename )  {
        int count = p_bins*v_bins;
        int * data = (int *) calloc(count,sizeof(int));
        double v_bin_size=(v_max-v_min)/v_bins;
        double p_bin_size=(p_max-p_min)/p_bins;
    for(CParIter pti(P,0);pti.isValid() ; ++pti){
        auto & particles = pti.GetArrayOfStructs();
        for (auto &p : particles){
            int vel_b = (int)((p.rdata(2+v_comp)-v_min)/v_bin_size);
            int pos_b = (int)((p.pos(p_comp)-p_min)/p_bin_size);    ;
            data[pos_b+vel_b*p_bins]+=1;
        }
    }
    amrex::ParallelDescriptor::ReduceIntSum(data,count);
    if (amrex::ParallelContext::IOProcessorSub()) {
    std::ofstream ofs(filename, std::ofstream::out);
    for(int j =0; j < v_bins; j++){
        for(int i =0;i< p_bins; i++){
            amrex::Print(ofs) << i*p_bin_size+p_min << " " << j* v_bin_size +v_min <<" " << data[i+j*p_bins] << "\n";
        }
    }
    ofs.close();
    }
    free(data);
}



#endif
