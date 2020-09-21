#ifndef amrex_util
#define amrex_util
#include "AMReX_Array.H"
#include "AMReX_Geometry.H"
#include "AMReX_MultiFab.H"
#include <AMReX_PlotFileUtil.H>
#include "particle_defs.hpp"
#include "w_defs.hpp"
#include <tuple>
#include <math.h>
// What cell index is a given point in?
// This is equivalent with the index for the "Lower left" corner
amrex::IntArray get_point_cell(const amrex::Geometry geom,const amrex::RealArray pos); 
void add_single_particle( CParticleTile&particlet ,amrex::RealArray pos , amrex::RealArray vel, double m,double q);
std::array<int,3> get_num_segments(const amrex::Geometry geom,const amrex::RealArray x_start,const amrex::RealArray  x_end);
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


template<int W_range>
void get_particle_number_density(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab &P_dens){
    P_dens.setVal(0); 
    P_dens.setBndry(0);
    for (amrex::MFIter mfi(P_dens); mfi.isValid(); ++mfi){
    
        // Each grid,tile has a their own local particle container
        auto& Part = P.GetParticles(0)[std::make_pair(mfi.index(),mfi.LocalTileIndex())];
        auto&  particles = Part.GetArrayOfStructs();
        amrex::Array4<amrex::Real> const& P_dens_loc = P_dens[mfi].array(); 

    const auto low = geom.ProbLo();
    const auto Ics = geom.InvCellSize();
    // Remote particle grid lower corner

     for(auto& p : particles){
        int coord[3];
        coord[X]=floor((p.pos(X) -low[X])*Ics[X]);
        coord[Y]=floor((p.pos(Y) -low[Y])*Ics[Y]);
        coord[Z]=floor((p.pos(Z) -low[Z])*Ics[Z]);
        
            

        auto px=(p.pos(X)-low[X])*Ics[X];
        auto py=(p.pos(Y)-low[Y])*Ics[Y];
        auto pz=(p.pos(Z)-low[Z])*Ics[Z];
        
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
                       P_dens_loc(nx,ny,nz,X)+=
                          Wp(px-nx)*Wp(py-ny)*Wp(pz-nz);
            }
        }
        }

    }

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
template<int comp>
int get_num_segments(const amrex::Geometry geom, amrex::Real x_start,amrex::Real x_end ){
    auto start_idx=get_point_line<comp>(geom,x_start);
    auto end_idx=get_point_line<comp>(geom,x_end);
    return abs(start_idx-end_idx)+1;

}


template<int comp> 
std::pair<amrex::Real,int> reflect_boundary(const amrex::Geometry geom,amrex::Real pos ){
    if(pos < geom.ProbLo(comp) +PART_BOUND*geom.CellSize(comp)){
            return {-pos + 2*(geom.ProbLo(comp) + PART_BOUND*geom.CellSize(comp) ),-1}; 
    }else if(pos > geom.ProbHi(comp) - PART_BOUND*geom.CellSize(comp) )
    {
            return {-pos + 2*(geom.ProbHi(comp) - PART_BOUND*geom.CellSize(comp) ),-1};
    }
    else{
        return {pos,1}; 
    }
}



// Only two segments
// system starts at (0,0) 
// cell size is 1 
template<int comp>
inline int construct_segments(amrex::Real x_start ,amrex::Real x_end,std::array<amrex::Real,3> const&seg_points,std::array<int,2> const& seg_idx){
    seg_idx[0]=floor(x_start);
    seg_idx[1]=floor(x_end);
    int diff=seg_idx[1]-seg_idx[0];
    int ng=abs(diff)+1;
    seg_points[0]=x_start;
    if(ng==1){
    seg_points[1]=seg_idx[0]+diff;
    seg_points[2]=x_end;
    }
    else{
    seg_points[1]=x_end; 
    }
    return ng; 
}




// List of segments with start,end , cell_idexx. all 1D
// How do we handle periodicity?
template<int comp>
void get_segment_list(const amrex::Geometry geom,
        std::array<std::tuple<amrex::Real,amrex::Real,int>,2> &segments,int num_segments, amrex::Real x_start ,amrex::Real x_end){
        
    const auto cellsize=geom.CellSize(comp);
    const auto problo = geom.ProbLo(comp);
    double segment_start = x_start;
    double segment_end;
    int segment_index=get_point_line<comp>(geom,segment_start);
        int sig;
    if(x_start > x_end ){
        sig=-1; 


    }else{
        sig=1;
    }
    if(num_segments == 2){
            segment_end = (segment_index+ (sig+1)/2)*cellsize + problo; 
            segments[0]=std::make_tuple(segment_start,segment_end,segment_index);
            segment_index =segment_index+sig;
        
            segment_start=(segment_index -((sig-1))/2)*cellsize+problo;
            if(!geom.isPeriodic(comp) && ((segment_index == geom.Domain().smallEnd(comp)+PART_BOUND ) || (segment_index == geom.Domain().bigEnd(comp)-PART_BOUND) )){
              auto res=reflect_boundary<comp>(geom,x_end);
              x_end=res.first;
              sig*=-1;
            }
        segment_end = x_end; 
        segments[1]=std::make_tuple(segment_start,segment_end,segment_index); 
    }
    else{
        segment_end = x_end; 
        segments[0]=std::make_tuple(segment_start,segment_end,segment_index);
    }
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
