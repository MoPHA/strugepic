#ifndef amrex_util
#define amrex_util
#include "AMReX_Array.H"
#include "AMReX_Geometry.H"
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
std::vector<std::tuple<amrex::Real,amrex::Real,int>> get_segment_list(const amrex::Geometry geom, amrex::Real x_start ,amrex::Real x_end){
        
    std::vector<std::tuple<amrex::Real,amrex::Real,int>> seg_list;
    const auto cellsize=geom.CellSize(comp);
    const auto problo = geom.ProbLo(comp);
    const auto probhi = geom.ProbHi(comp);
    const auto num_segments = get_num_segments<comp>(geom,x_start,x_end);
    double segment_start = x_start;
    double segment_end;
    int segment_index=get_point_line<comp>(geom,segment_start);
        int sig;
    if(x_start > x_end ){
        sig=-1; 


    }else{
        sig=1;
    }
    for(int seg=1;seg < num_segments;seg++){ 
            segment_end = (segment_index+ (sig+1)/2)*cellsize + problo; 
            seg_list.push_back(std::make_tuple(segment_start,segment_end,segment_index));
            segment_index =segment_index+sig;//(((segment_index +sig) % index_mod) +index_mod) %index_mod ;
        
            segment_start=(segment_index -((sig-1))/2)*cellsize+problo;
            if(!geom.isPeriodic(comp) && ((segment_index == geom.Domain().smallEnd(comp)+PART_BOUND ) || (segment_index == geom.Domain().bigEnd(comp)-PART_BOUND) )){
              auto res=reflect_boundary<comp>(geom,x_end);
              x_end=res.first;
              sig*=-1;
            }

    }

   segment_end = x_end;// wrapMinMax(x_end,problo,probhi);
   
   seg_list.push_back(std::make_tuple(segment_start,segment_end,segment_index));
   return seg_list;

}

std::pair<amrex::Real,amrex::Real> get_total_energy(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab &E, amrex::MultiFab &B );
std::pair<std::array<amrex::Real,3>,std::array<amrex::Real,3>> get_total_momentum(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab &E, amrex::MultiFab &B);
void add_particle_n_per_cell(const amrex::Geometry geom, CParticleContainer&P,double m,double q,double v,int n);


#endif
