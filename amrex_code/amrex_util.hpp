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

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}




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


// List of segments with start,end , cell_idexx. all 1D
// How do we handle periodicity?
template<int comp>
std::vector<std::tuple<amrex::Real,amrex::Real,int>> get_segment_list(const amrex::Geometry geom, amrex::Real x_start ,amrex::Real x_end){
        
    std::vector<std::tuple<amrex::Real,amrex::Real,int>> seg_list;
    const auto cellsize=geom.CellSize(comp);
    const auto problo = geom.ProbLo(comp);
    // Periodic index
    int index_mod = geom.Domain().hiVect()[comp]+1;
    const auto prob_size=geom.ProbHi(comp)-problo;
    auto num_segments = get_num_segments<comp>(geom,x_start,x_end);
    

    double segment_start = x_start;
    double segment_end;
    int segment_index=get_point_line<comp>(geom,segment_start);
    
    if(x_start > x_end ){
            for(int seg=1;seg < num_segments;seg++){ 
                segment_end = (segment_index)*cellsize + problo; 
                seg_list.push_back(std::make_tuple(segment_start,segment_end,segment_index));
                segment_index =(segment_index -1) % index_mod ;
                segment_start=fmod(segment_end-problo,prob_size)+problo;

            }

        }else{
            for(int seg=1;seg < num_segments;seg++){ 
                segment_end = (segment_index+1)*cellsize + problo; 
                seg_list.push_back(std::make_tuple(segment_start,segment_end,segment_index));
                segment_index =(segment_index +1)% index_mod;
                segment_start=fmod(segment_end-problo,prob_size)+problo;

            }
        }
        segment_end = fmod(x_end-problo,prob_size)+problo;;
        seg_list.push_back(std::make_tuple(segment_start,segment_end,segment_index));
        return seg_list;

}

#endif
