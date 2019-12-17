#ifndef PROPAGATOR
#define PROPAGATOR
#include<AMReX_Geometry.H>
#include "../include/interpolation.hpp"
#include "particle_defs.hpp"
#include "amrex_util.hpp"

// These are All local update functions

void push_B( amrex::Box const& bx,  amrex::Array4<amrex::Real> const& B, amrex::Array4<amrex::Real> const& E,double dt);

void push_E( amrex::Box const& bx,  amrex::Array4<amrex::Real> const& E, amrex::Array4<amrex::Real> const& B,double dt);

void push_V_E( CParticles&particles, const amrex::Geometry geom,amrex::Array4<amrex::Real> const& E ,double dt);

void Theta_E(const amrex::Geometry geom,amrex::Box const& bx,amrex::Array4<amrex::Real> const& E,amrex::Array4<amrex::Real> const& B,CParticles&particles,double dt );

void Theta_B(amrex::Box const& bx,amrex::Array4<amrex::Real> const& E,amrex::Array4<amrex::Real> const& B,double dt );

void Theta_x(CParticleContainer&ParticleC, CParticles&local_particles,const CNParticles&neighbour_particles, const amrex::Geometry geom,amrex::Array4<amrex::Real> const& E, amrex::Array4<amrex::Real> const& B ,double dt);

void Theta_y(CParticleContainer&ParticleC, CParticles&local_particles,const CNParticles&neighbour_particles, const amrex::Geometry geom,amrex::Array4<amrex::Real> const& E, amrex::Array4<amrex::Real> const& B ,double dt);

void Theta_z(CParticleContainer&ParticleC, CParticles&local_particles,const CNParticles&neighbour_particles, const amrex::Geometry geom,amrex::Array4<amrex::Real> const& E, amrex::Array4<amrex::Real> const& B ,double dt);

// Global update, they also handle triggering the updates
void G_Theta_E(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab &E, amrex::MultiFab &B,double dt );
void G_Theta_B(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab &E, amrex::MultiFab &B,double dt );
void G_Theta_x(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab &E, amrex::MultiFab &B,double dt );
void G_Theta_y(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab &E, amrex::MultiFab &B,double dt );
void G_Theta_z(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab &E, amrex::MultiFab &B,double dt );





// This is for one particle type
// If there are several you need to do this again
// 0 -> x  , 1-> y 2->z
// Local and neighbour particle list are not the same data structure
template<int comp,class Pcontainter>
void push_E_part(const Pcontainter&particles, const amrex::Geometry geom,amrex::Array4<amrex::Real> const& E ,double dt){ 
    const int idx_list[4]={-1,0,1,2};
    const  double coef = particles[0].rdata(1)*geom.InvCellSize(0)*geom.InvCellSize(1)*geom.InvCellSize(2);
    const auto low = geom.ProbLo();
    const auto Ics = geom.InvCellSize();

    for(auto& p : particles){
        const auto p_segments = get_segment_list<comp>(geom, p.pos(comp) , p.pos(comp)+dt*p.rdata(comp+2));
        auto coord =get_point_cell(geom,{p.pos(0),p.pos(1),p.pos(2)}) ;
        for(auto seg: p_segments){
            coord[comp] = std::get<2>(seg);
            auto i_s = std::get<0>(seg); 
            auto i_e = std::get<1>(seg);

    for(auto k: idx_list){
            for(auto j: idx_list){
                for(auto i: idx_list){
                    amrex::Real res;
                    auto cx = coord[0]+i;
                    auto cy = coord[1]+j;
                    auto cz = coord[2]+k;
                    using namespace strugepic;

                    if(E.contains(cx,cy,cz)){
                    switch(comp){
                        case 0:
                            res=I_W12((i_s-low[0])*Ics[0]- cx ,(i_e-low[0])*Ics[0]- cx)*W1((p.pos(1)-low[1])*Ics[1]-cy  )*W1((p.pos(2)-low[2])*Ics[2]-cz );
                            break;
                        case 1:
                            res=I_W12((i_s-low[1])*Ics[1]- cy ,(i_e-low[1])*Ics[1]- cy)*W1((p.pos(2)-low[2])*Ics[2]-cz )*W1( (p.pos(0)-low[0])*Ics[0]-cx   );
                            break;
                        case 2:
                            res=I_W12((i_s-low[2])*Ics[2]- cz ,(i_e-low[2])*Ics[2]- cz)*W1( (p.pos(0)-low[0])*Ics[0]-cx   )*W1((p.pos(1)-low[1])*Ics[1]-cy  );
                            break;
                
                    }
                    E(cx,cy,cz,comp)-=coef*p.rdata(comp+2)*res;
                    }
                    
                }

            }
    }
    }
}
    
}


template<int comp>
void push_E_p(const CParticles&local_particles,const CNParticles&neighbour_particles, const amrex::Geometry geom,amrex::Array4<amrex::Real> const& E ,double dt){
    push_E_part<comp,CParticles>(local_particles,geom,E,dt); 
    if(neighbour_particles.size()!=0){
    push_E_part<comp,CNParticles>(neighbour_particles,geom,E,dt); 
    }
}

template<int comp>
void push_V_p(CParticleContainer&p_container ,CParticles&local_particles,double dt){
        for(auto &p: local_particles){
            p.pos(comp)+=dt*p.rdata(2+comp);
            p_container.Reset(p,true);
        }
}


template<int comp>
void push_B_p(CParticles&particles, const amrex::Geometry geom, const amrex::Array4<amrex::Real> & B ,double dt){ 
    const int idx_list[4]={-1,0,1,2};
    const  double coef = particles[0].rdata(1)/particles[0].rdata(0);
    const auto low = geom.ProbLo();
    const auto Ics = geom.InvCellSize();

    for(auto& p : particles){
        const auto p_segments = get_segment_list<comp>(geom, p.pos(comp) , p.pos(comp)+dt*p.rdata(comp+2));
        auto coord =get_point_cell(geom,{p.pos(0),p.pos(1),p.pos(2)}) ;

        amrex::Real res_c1=0;
        amrex::Real res_c2=0;
        for(auto seg: p_segments){
            coord[comp] = std::get<2>(seg);
            auto i_s = std::get<0>(seg); 
            auto i_e = std::get<1>(seg);



// The W_1^{(2)} function is supported over -1 < x< 2 so symmetric looping is not the best, but will do for now.
    double norm_res=0;
    for(auto k: idx_list){
            for(auto j: idx_list){
                for(auto i: idx_list){
                    auto cx = coord[0]+i;
                    auto cy = coord[1]+j;
                    auto cz = coord[2]+k;
                    using namespace strugepic;
 //                   norm_res+=W( (p.pos(0)-low[0])*Ics[0]-cx  , (p.pos(1)-low[1])*Ics[1]-cy    ,  (p.pos(2)-low[2])*Ics[2]-cz ) ;
                    switch(comp){
                        case 0:
                            res_c1+=B(cx,cy,cz,1)*I_W12((i_s-low[0])*Ics[0]- cx ,(i_e-low[0])*Ics[0]- cx)*W1(  (p.pos(1)-low[1])*Ics[1]-cy  )*W12( (p.pos(2)-low[2])*Ics[2]-cz  );
                            res_c2-=B(cx,cy,cz,2)*I_W12((i_s-low[0])*Ics[0]- cx ,(i_e-low[0])*Ics[0]- cx)*W12( (p.pos(1)-low[1])*Ics[1]-cy  )*W1( (p.pos(2)-low[2])*Ics[2]-cz  );
                            break;
                        case 1:
                            res_c1+=B(cx,cy,cz,2)*I_W12((i_s-low[1])*Ics[1]- cy ,(i_e-low[1])*Ics[1]- cy)*W12( (p.pos(0)-low[0])*Ics[0]-cx  )*W1( (p.pos(2)-low[2])*Ics[2]-cz  );
                            res_c2-=B(cx,cy,cz,0)*I_W12((i_s-low[1])*Ics[1]- cy ,(i_e-low[1])*Ics[1]- cy)*W1( (p.pos(0)-low[0])*Ics[0]-cx  )*W12( (p.pos(2)-low[2])*Ics[2]-cz  );
                            break;
                        case 2:
                            res_c1+=B(cx,cy,cz,0)*I_W12((i_s-low[2])*Ics[2]- cz ,(i_e-low[2])*Ics[2]- cz)*W1( (p.pos(0)-low[0])*Ics[0]-cx  )*W12((p.pos(1)-low[1])*Ics[1]-cy  );
                            res_c2-=B(cx,cy,cz,1)*I_W12((i_s-low[2])*Ics[2]- cz ,(i_e-low[2])*Ics[2]- cz)*W12( (p.pos(0)-low[0])*Ics[0]-cx  )*W1((p.pos(1)-low[1])*Ics[1]-cy);
                            break;
                
                    }

                    
                }

            }
    }   
    }
        p.rdata( (comp +2 )% 3 +2   )+=coef*res_c1*p.rdata(comp+2); 
        p.rdata( (comp+1) % 3 +2  )+=coef*res_c2*p.rdata(comp+2); 
}
    
}


#endif

