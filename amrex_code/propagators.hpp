#ifndef PROPAGATOR
#define PROPAGATOR
#include<AMReX_Geometry.H>
#include "../include/interpolation.hpp"
#include "particle_defs.hpp"
#include "amrex_util.hpp"

// These are All local update functions, I.e they operate only on local data 

void push_B_E(const amrex::Geometry geom, amrex::Box const& bx,  amrex::Array4<amrex::Real> const& B, amrex::Array4<amrex::Real> const& E,double dt);

void push_E_B(const amrex::Geometry geom, amrex::Box const& bx,  amrex::Array4<amrex::Real> const& E, amrex::Array4<amrex::Real> const& B,double dt);

void push_V_E( CParticles&particles, const amrex::Geometry geom,amrex::Array4<amrex::Real> const& E ,double dt);

void Theta_E(const amrex::Geometry geom,amrex::Box const& bx,amrex::Array4<amrex::Real> const& E,amrex::Array4<amrex::Real> const& B,CParticles&particles,double dt );

void Theta_B(const amrex::Geometry geom,amrex::Box const& bx,amrex::Array4<amrex::Real> const& E,amrex::Array4<amrex::Real> const& B,double dt );


// Global update, they also handle triggering global communication
void G_Theta_E(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab &E, amrex::MultiFab &B,double dt );
void G_Theta_B(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab &E, amrex::MultiFab &B,double dt );






// This is for one particle type
// If there are several you need to do this again
// 0 -> x  , 1-> y 2->z
// Local and neighbour particle list are not the same data structure
template<int comp,class Pcontainter>
void push_E_part(const Pcontainter&particles, const amrex::Geometry geom,amrex::Array4<amrex::Real> const& E ,double dt){ 
    const int idx_list[4]={-1,0,1,2};
    const  double coef = particles[0].rdata(Q)*geom.InvCellSize(X)*geom.InvCellSize(Y)*geom.InvCellSize(Z);
    const auto low = geom.ProbLo();
    const auto Ics = geom.InvCellSize();

    for(auto& p : particles){
        const auto p_segments = get_segment_list<comp>(geom, p.pos(comp) , p.pos(comp)+dt*p.rdata(comp+2));
        auto coord =get_point_cell(geom,{p.pos(X),p.pos(Y),p.pos(Z)}) ;
        for(auto seg: p_segments){
            coord[comp] = std::get<2>(seg);
            auto i_s = std::get<0>(seg); 
            auto i_e = std::get<1>(seg);

            for(auto k: idx_list){
                for(auto j: idx_list){
                    for(auto i: idx_list){
                        amrex::Real res;
                        int cl[3];
                        cl[X] = coord[X]+i;
                        cl[Y] = coord[Y]+j;
                        cl[Z] = coord[Z]+k;
                        auto comp_u = (comp+1)%3;
                        auto comp_l = (comp+2)%3;
                        
                        using namespace strugepic;
                        if(E.contains(cl[X],cl[Y],cl[Z])){
                            res=I_W12((i_s-low[comp])*Ics[comp]- cl[comp] ,(i_e-low[comp])*Ics[comp]- cl[comp])
                                *W1((p.pos(comp_u)-low[comp_u])*Ics[comp_u]-cl[comp_u])
                                *W1((p.pos(comp_l)-low[comp_l])*Ics[comp_l]-cl[comp_l]);

                            E(cl[X],cl[Y],cl[Z],comp)-=coef*res;
                        }
                    }
                }
            }
        }
    }
}


template<int comp>
void push_E_pos(const CParticles&local_particles,const CNParticles&neighbour_particles, const amrex::Geometry geom,amrex::Array4<amrex::Real> const& E ,double dt){
    push_E_part<comp,CParticles>(local_particles,geom,E,dt); 
    if(neighbour_particles.size()!=0){
    push_E_part<comp,CNParticles>(neighbour_particles,geom,E,dt); 
    }
}

template<int comp>
void push_pos_pos(CParticleContainer&p_container ,CParticles&local_particles,double dt){
        for(auto &p: local_particles){
            p.pos(comp)+=dt*p.rdata(2+comp);
            p_container.Reset(p,true);
        }
}


template<int comp>
void push_B_pos(CParticles&particles, const amrex::Geometry geom, const amrex::Array4<amrex::Real> & B ,double dt){ 
    const int idx_list[4]={-1,0,1,2};
    const auto low = geom.ProbLo();
    const auto Ics = geom.InvCellSize();
    // A scaling factor is needed here if dx,dy,dz are not equal!
    // We choose to scale to units of dx
    const  double coef = particles[0].rdata(Q)/particles[0].rdata(M)*geom.CellSize(comp);

    for(auto& p : particles){
        const auto p_segments = get_segment_list<comp>(geom, p.pos(comp) , p.pos(comp)+dt*p.rdata(comp+2));
        auto coord =get_point_cell(geom,{p.pos(X),p.pos(Y),p.pos(Z)}) ;
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
                        int cl[3];
                        cl[X] = coord[X]+i;
                        cl[Y] = coord[Y]+j;
                        cl[Z] = coord[Z]+k;
                        auto comp_u = (comp+1)%3;
                        auto comp_l = (comp+2)%3;
                        
                        using namespace strugepic;//norm_res+=W( (p.pos(0)-low[0])*Ics[0]-cx  , (p.pos(1)-low[1])*Ics[1]-cy    ,  (p.pos(2)-low[2])*Ics[2]-cz ) ;            
                        //const double B_pseudo[3]={0,0,0.05};

                        //res_c1+=B_pseudo[comp_u]
                        res_c1+=B(cl[X],cl[Y],cl[Z],comp_u)
                        *I_W12((i_s-low[comp])*Ics[comp]- cl[comp] ,(i_e-low[comp])*Ics[comp]- cl[comp])
                        *W1(  (p.pos(comp_u)-low[comp_u])*Ics[comp_u]-cl[comp_u])
                        *W12( (p.pos(comp_l)-low[comp_l])*Ics[comp_l]-cl[comp_l]);
                        
                        //res_c2-=B_pseudo[comp_l]
                        res_c2-=B(cl[X],cl[Y],cl[Z],comp_l)
                        *I_W12((i_s-low[comp])*Ics[comp]- cl[comp] ,(i_e-low[comp])*Ics[comp]- cl[comp])
                        *W12( (p.pos(comp_u)-low[comp_u])*Ics[comp_u]-cl[comp_u])
                        *W1(  (p.pos(comp_l)-low[comp_l])*Ics[comp_l]-cl[comp_l]);

                    } // over i
                } // over j
            } // over k  
        } // over segments

        p.rdata( (comp +2 )% 3 +2   )+=coef*res_c1;
        p.rdata( (comp+1) % 3 +2  )+=coef*res_c2; 

//        p.rdata(VX) = sgn(p.rdata(VX))*std::min(fabs(p.rdata(VX)),0.9);
//        p.rdata(VY) = sgn(p.rdata(VY))*std::min(fabs(p.rdata(VY)),0.9);
//        p.rdata(VZ) = sgn(p.rdata(VZ))*std::min(fabs(p.rdata(VZ)),0.9);
    } // over particles
} // function end


template <int coord>
void Theta(CParticleContainer&ParticleC, CParticles&local_particles,const CNParticles&neighbour_particles, const amrex::Geometry geom,amrex::Array4<amrex::Real> const& E, amrex::Array4<amrex::Real> const& B ,double dt)
{
    push_E_pos<coord>(local_particles,neighbour_particles,geom,E,dt);
    push_B_pos<coord>(local_particles,geom,B,dt);
    push_pos_pos<coord>(ParticleC,local_particles,dt);
}



template <int coord>
void G_Theta(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab &E, amrex::MultiFab &B,double dt ){

    for (CParIter pti(P, 0); pti.isValid(); ++pti) {
        auto&  particles = pti.GetArrayOfStructs();
        const auto& n_particles = P.GetNeighbors(0,pti.index(),pti.LocalTileIndex());
        amrex::FArrayBox& efab = E[pti];
        amrex::FArrayBox& bfab =B[pti];
        amrex::Array4<amrex::Real> const& E_loc = efab.array();
        amrex::Array4<amrex::Real> const& B_loc = bfab.array(); 

        Theta<coord>(P,particles,n_particles,geom,E_loc,B_loc,dt);
    }
    P.Redistribute();
    P.updateNeighbors();
    E.FillBoundary(geom.periodicity());
    B.FillBoundary(geom.periodicity());
    
}

#endif

