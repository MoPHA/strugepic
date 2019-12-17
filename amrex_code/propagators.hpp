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
    const  double dx = geom.CellSize(0);
    const  double dy = geom.CellSize(1);
    const  double dz = geom.CellSize(2);

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

                    if(E.contains(coord[0]+i,coord[1]+j,coord[2]+k)){
                    switch(comp){
                        case 0:
                            res=strugepic::I_W12(i_s - (coord[0]+i)*dx,i_e - (coord[0]+i)*dx)*strugepic::W1(p.pos(1) - (coord[1]+j)*dy  )*strugepic::W1(p.pos(2) - (coord[2]+k)*dz  );
                            break;
                        case 1:
                            res=strugepic::I_W12(i_s- (coord[1]+j)*dy,i_e - (coord[1]+j)*dy)*strugepic::W1(p.pos(2)- (coord[2]+k)*dz  )*strugepic::W1(p.pos(0) - (coord[0]+i)*dx  );
                            break;
                        case 2:
                            res=strugepic::I_W12(i_s- (coord[2]+k)*dz,i_e - (coord[2]+k)*dz)*strugepic::W1(p.pos(0)- (coord[0]+i)*dx  )*strugepic::W1(p.pos(1) - (coord[1]+j)*dy  );
                            break;
                
                    }
                    E(coord[0]+i,coord[1]+j,coord[2]+k,comp)-=coef*p.rdata(comp+2)*res;
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
    const  double dx = geom.CellSize(0);
    const  double dy = geom.CellSize(1);
    const  double dz = geom.CellSize(2);

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
    for(auto k: idx_list){
            for(auto j: idx_list){
                for(auto i: idx_list){
                    switch(comp){
                        case 0:
                            res_c1+=B(coord[0]+i,coord[1]+j,coord[2]+k,1)*strugepic::I_W12(i_s- (coord[0]+i)*dx,i_e - (coord[0]+i)*dx)*strugepic::W1(p.pos(1)- (coord[1]+j)*dy  )*strugepic::W12(p.pos(2) - (coord[2]+k)*dz  );
                            res_c2-=B(coord[0]+i,coord[1]+j,coord[2]+k,2)*strugepic::I_W12(i_s- (coord[0]+i)*dx,i_e - (coord[0]+i)*dx)*strugepic::W12(p.pos(1)- (coord[1]+j)*dy  )*strugepic::W1(p.pos(2) - (coord[2]+k)*dz  );
                            break;
                        case 1:
                            res_c1+=B(coord[0]+i,coord[1]+j,coord[2]+k,2)*strugepic::I_W12(i_s- (coord[1]+j)*dy,i_e - (coord[1]+j)*dy)*strugepic::W12(p.pos(0)- (coord[0]+i)*dx  )*strugepic::W1(p.pos(2) - (coord[2]+k)*dz  );
                            res_c2-=B(coord[0]+i,coord[1]+j,coord[2]+k,0)*strugepic::I_W12(i_s - (coord[1]+j)*dy,i_e - (coord[1]+j)*dy)*strugepic::W1(p.pos(0) - (coord[0]+i)*dx  )*strugepic::W12(p.pos(2) - (coord[2]+k)*dz  );
                            break;
                        case 2:
                            res_c1+=B(coord[0]+i,coord[1]+j,coord[2]+k,0)*strugepic::I_W12(i_s - (coord[2]+k)*dz,i_e - (coord[2]+k)*dz)*strugepic::W1(p.pos(0) - (coord[0]+i)*dx  )*strugepic::W12(p.pos(1) - (coord[1]+j)*dy  );
                            res_c2-=B(coord[0]+i,coord[1]+j,coord[2]+k,1)*strugepic::I_W12(i_s- (coord[2]+k)*dz,i_e - (coord[2]+k)*dz)*strugepic::W12(p.pos(0)- (coord[0]+i)*dx  )*strugepic::W1(p.pos(1) - (coord[1]+j)*dy  );
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

