#ifndef amrex_util
#define amrex_util
#include "AMReX_Array.H"
#include "AMReX_Geometry.H"
#include "AMReX_MultiFab.H"
#include "particle_defs.hpp"
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
void get_particle_number_density(const amrex::Geometry geom,const amrex::Geometry aux_geom,CParticleContainer&P, amrex::MultiFab &P_dens,amrex::MultiFab & P_dens_aux);
double uniform_density(const amrex::Geometry geom,int i ,int j ,int k);
double gaussian_dist(double pos,double simu_len);
void set_field_gradient_gaussian_x(amrex::MultiFab &A,double A_max,double simu_len);

template<int comp>
void shift_and_grow(amrex::BoxArray &ba, int nghost){
    // Create a list of 
    // box index and coordinate of lower left corner
    // no-overlaps in the input box array
    std::vector<std::pair<int,int>> cc;
    for(int i=0;i<ba.size();i++){
        cc.push_back({ba[i].loVect3d()[comp],i});

    }
    std::sort(cc.begin(),cc.end());
    
    // Shift the coordinates
    int num_box=0;
    int cval=cc[0].first;
    for(auto e:cc){
       if(e.first > cval){
            num_box+=1;
            cval=e.first;
       } 
        ba.set(e.second,ba[e.second].shift(comp,2*nghost*num_box).grow(comp,nghost));

    }
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
    SimulationIO(amrex::Geometry geom,amrex::Geometry aux_geom,amrex::BoxArray gba,amrex::MultiFab & E,amrex::MultiFab & B,CParticleContainer &P,double dt,std::string data_folder_name);
    void write(int step,bool checkpoint=false,bool particles=false);
    void read(int step);
    void dump_pdens(std::string filename);
    void dump_E_field(std::string filename);
    void dump_B_field(std::string filename);
    private:
        amrex::Geometry geom;
        amrex::Geometry aux_geom;
        amrex::MultiFab & E;
        amrex::MultiFab & B;
        amrex::MultiFab  Pdens;
        amrex::MultiFab  Pdens_aux;
        CParticleContainer &P;
        double dt;
        std::string data_folder_name; 


};



double wrapMax(double x, double max);
double wrapMinMax(double x , double min ,double max);
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
            segment_index =segment_index+sig;//(((segment_index +sig) % index_mod) +index_mod) %index_mod ;
        
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
std::pair<std::array<amrex::Real,3>,std::array<amrex::Real,3>> get_total_momentum(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab &E, amrex::MultiFab &B);
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
void set_temprature_gradient(CParticleContainer &P,double (*grad_fun) (double,double) ,double v_min,double v_max,double system_len){
    std::random_device rd;
    std::mt19937 mt(rd());
    
    double v_min_var=v_min*v_min;
    double v_max_var=v_max*v_max;
    double v_diff_var=v_max_var-v_min_var;
    std::uniform_real_distribution<double> angle(0,2*M_PI);
    for(CParIter pti(P,0);pti.isValid() ; ++pti){ 
        auto & particles = pti.GetArrayOfStructs();
        for (auto &p : particles){
            double v = sqrt(v_min_var+grad_fun(p.pos(comp),system_len)*v_diff_var);
            double ang = angle(mt);
            std::normal_distribution<double> vel(0,v);
            double vx = vel(mt);
            double vy = vel(mt);
            double v_mag=sqrt(vx*vx+vy*vy);
            double rad = p.rdata(M)*v_mag/(p.rdata(Q)*2678.2);

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
