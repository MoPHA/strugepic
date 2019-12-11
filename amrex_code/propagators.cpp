#include "AMReX_Geometry.H"
#include "particle_defs.hpp"
#include "propagators.hpp"
#include "amrex_util.hpp"


void Theta_E(const amrex::Geometry geom,amrex::Box const& bx,amrex::Array4<amrex::Real> const& E,amrex::Array4<amrex::Real> const& B,CParticles&particles,double dt ){
    push_B(bx, B,E,dt);
    push_V_E(particles,geom,E,dt);
}

void Theta_B(amrex::Box const& bx,amrex::Array4<amrex::Real> const& E,amrex::Array4<amrex::Real> const& B,double dt ){
    push_E(bx,E,B,dt);
}



void push_B( amrex::Box const& bx,  amrex::Array4<amrex::Real> const& B, amrex::Array4<amrex::Real> const& E,double dt){
   const auto lo = amrex::lbound(bx);
   const auto hi = amrex::ubound(bx);
   for     (int k = lo.z; k <= hi.z; ++k) {
     for   (int j = lo.y; j <= hi.y; ++j) {
       for (int i = lo.x; i <= hi.x; ++i) { 
         B(i,j,k,0) -=dt* (E(i+1,j,k,0)-E(i,j,k,0));
         B(i,j,k,1) -=dt* (E(i,j+1,k,1)-E(i,j,k,1));
         B(i,j,k,2) -=dt* (E(i,j,k+1,2)-E(i,j,k,2));
       }
     }
   }

}

void push_E( amrex::Box const& bx,  amrex::Array4<amrex::Real> const& E, amrex::Array4<amrex::Real> const& B,double dt){
   const auto lo = amrex::lbound(bx);
   const auto hi = amrex::ubound(bx);
   for     (int k = lo.z; k <= hi.z; ++k) {
     for   (int j = lo.y; j <= hi.y; ++j) {
       for (int i = lo.x; i <= hi.x; ++i) { 
         E(i,j,k,0) +=dt* (B(i+1,j,k,0)-B(i,j,k,0));
         E(i,j,k,1) +=dt* (B(i,j+1,k,1)-B(i,j,k,1));
         E(i,j,k,2) +=dt* (B(i,j,k+1,2)-B(i,j,k,2));
       }
     }
   }

}


void push_V_E( CParticles&particles, const amrex::Geometry geom,amrex::Array4<amrex::Real> const& E ,double dt){

    // Basis functions are supported over two cells in each direction
    const int idx_list[4]={-1,0,1,2};
    // All the particles in a give container have the same mass and charge
    // We could probably do without saving the mass and charge for each particle
    const double m= particles[0].rdata(0);
    const double q= particles[0].rdata(1);
    const double coef = dt*q/m;
    const double inv_dx = geom.InvCellSize(0);
    const double inv_dy = geom.InvCellSize(1);
    const double inv_dz = geom.InvCellSize(2);
    
    for(auto& p : particles){
        auto index =get_point_cell(geom,{p.pos(0),p.pos(1),p.pos(2)}) ;
       
        double dvx=0;
        double dvy=0;
        double dvz=0;


        for(auto k: idx_list){
            for(auto j: idx_list){
                for(auto i: idx_list){
                  auto W = strugepic::W_s1({ p.pos(0)*inv_dx - (index[0]+i),p.pos(1)*inv_dy -(index[1]+j),p.pos(2)*inv_dz -(index[2]+k) }); 
                   dvx+= E(index[0]+i,index[1]+j,index[2]+k,0)*W[0]; 
                   dvy+= E(index[0]+i,index[1]+j,index[2]+k,1)*W[1];
                   dvz+= E(index[0]+i,index[1]+j,index[2]+k,2)*W[2];
                }
            }
        }
        p.rdata(2)+=dvx*coef;
        p.rdata(3)+=dvy*coef;
        p.rdata(4)+=dvz*coef;
    }
}

// This is for one particle type
// If there are several you need to do this again
// 0 -> x  , 1-> y 2->z
template<int comp>
void push_E(CParticles&particles, const amrex::Geometry geom,amrex::Array4<amrex::Real> const& E ,double dt){
    const  double coef = particles[0].rdata(1)*geom.InvCellSize(0)*geom.InvCellSize(1)*geom.InvCellSize(2);
    double dE =0;

    for(auto& p : particles){
    }

    dE*coef;
};
