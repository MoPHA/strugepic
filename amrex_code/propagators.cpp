#include "AMReX_Array.H"
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

void Theta_x(CParticleContainer&ParticleC, CParticles&local_particles,const CNParticles&neighbour_particles, const amrex::Geometry geom,amrex::Array4<amrex::Real> const& E, amrex::Array4<amrex::Real> const& B ,double dt){
    push_E_p<0>(local_particles,neighbour_particles,geom,E,dt);
    push_B_p<0>(local_particles,geom,B,dt);
    push_V_p<0>(ParticleC,local_particles,dt);
}
void Theta_y(CParticleContainer&ParticleC, CParticles&local_particles,const CNParticles&neighbour_particles, const amrex::Geometry geom,amrex::Array4<amrex::Real> const& E, amrex::Array4<amrex::Real> const& B ,double dt){
    push_E_p<1>(local_particles,neighbour_particles,geom,E,dt);
    push_B_p<1>(local_particles,geom,B,dt);
    push_V_p<1>(ParticleC,local_particles,dt);
}
void Theta_z(CParticleContainer&ParticleC, CParticles&local_particles,const CNParticles&neighbour_particles, const amrex::Geometry geom,amrex::Array4<amrex::Real> const& E, amrex::Array4<amrex::Real> const& B ,double dt){
    push_E_p<2>(local_particles,neighbour_particles,geom,E,dt);
    push_B_p<2>(local_particles,geom,B,dt);
    push_V_p<2>(ParticleC,local_particles,dt);
}


void push_B( amrex::Box const& bx,  amrex::Array4<amrex::Real> const& B, amrex::Array4<amrex::Real> const& E,double dt){
   const auto lo = amrex::lbound(bx);
   const auto hi = amrex::ubound(bx);
   for     (int k = lo.z; k <= hi.z; ++k) {
     for   (int j = lo.y; j <= hi.y; ++j) {
       for (int i = lo.x; i <= hi.x; ++i) { 
         B(i,j,k,0) -=dt*((E(i,j+1,k,2)-E(i,j,k,2))-( E(i,j,k+1,1)-E(i,j,k,1)));  
         B(i,j,k,1) -=dt*((E(i,j,k+1,0)-E(i,j,k,0))-( E(i+1,j,k,2)-E(i,j,k,2)));
         B(i,j,k,2) -=dt*((E(i+1,j,k,1)-E(i,j,k,1))-( E(i,j+1,k,0)-E(i,j,k,0)));
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
         E(i,j,k,0) +=dt*((B(i,j+1,k,2)-B(i,j,k,2))-( B(i,j,k+1,1)-B(i,j,k,1)));  
         E(i,j,k,1) +=dt*((B(i,j,k+1,0)-B(i,j,k,0))-( B(i+1,j,k,2)-B(i,j,k,2)));
         E(i,j,k,2) +=dt*((B(i+1,j,k,1)-B(i,j,k,1))-( B(i,j+1,k,0)-B(i,j,k,0)));
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
    const double dx = geom.CellSize(0);
    const double dy = geom.CellSize(1);
    const double dz = geom.CellSize(2);
    
    for(auto& p : particles){
        auto index =get_point_cell(geom,{p.pos(0),p.pos(1),p.pos(2)}) ;
       
        double dvx=0;
        double dvy=0;
        double dvz=0;


        for(auto k: idx_list){
            for(auto j: idx_list){
                for(auto i: idx_list){
                  auto W = strugepic::W_s1({ p.pos(0) - (index[0]+i)*dx,p.pos(1) -(index[1]+j)*dy,p.pos(2) -(index[2]+k)*dz }); 
                   dvx+= E(index[0]+i,index[1]+j,index[2]+k,0)*W[0]; 
                   dvy+= E(index[0]+i,index[1]+j,index[2]+k,1)*W[1];
                   dvz+= E(index[0]+i,index[1]+j,index[2]+k,2)*W[2];
                }
            }
        }
       // std::cout << "Speed updates from E: "<<dvx*coef  << "," << dvy*coef  << ","<< dvx*coef <<std::endl;
        p.rdata(2)+=dvx*coef;
        p.rdata(3)+=dvy*coef;
        p.rdata(4)+=dvz*coef;
    }
}


void G_Theta_E(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab&E, amrex::MultiFab&B,double dt ){


    for (CParIter pti(P, 0); pti.isValid(); ++pti) {
        auto&  particles = pti.GetArrayOfStructs();
        amrex::FArrayBox& efab = E[pti];
        amrex::Array4<amrex::Real> const& E_loc = efab.array();
         
        push_V_E(particles,geom,E_loc,dt);
    }

    for (amrex::MFIter mfi(E); mfi.isValid(); ++mfi){
        const amrex::Box& box = mfi.validbox();
        amrex::FArrayBox& fab = E[mfi];
        amrex::FArrayBox& fabB = B[mfi];
        amrex::Array4<amrex::Real> const& E_loc = fab.array();
        amrex::Array4<amrex::Real> const& B_loc = fabB.array(); 
        push_B(box, B_loc,E_loc,dt);

    }

    P.updateNeighbors();
    B.FillBoundary(geom.periodicity());


}



void G_Theta_B(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab&E, amrex::MultiFab&B,double dt ){

    for (amrex::MFIter mfi(E); mfi.isValid(); ++mfi){
        const amrex::Box& box = mfi.validbox();
        amrex::FArrayBox& fab = E[mfi];
        amrex::FArrayBox& fabB = B[mfi];
        amrex::Array4<amrex::Real> const& E_loc = fab.array();
        amrex::Array4<amrex::Real> const& B_loc = fabB.array();
    
        Theta_B(box,E_loc,B_loc,dt);

    }

    P.updateNeighbors();
    E.FillBoundary();

}
void G_Theta_x(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab&E, amrex::MultiFab&B,double dt ){
    for (CParIter pti(P, 0); pti.isValid(); ++pti) {
        auto&  particles = pti.GetArrayOfStructs();
        const auto& n_particles = P.GetNeighbors(0,pti.index(),pti.LocalTileIndex());
        amrex::FArrayBox& efab = E[pti];
        amrex::FArrayBox& bfab =B[pti];
        amrex::Array4<amrex::Real> const& E_loc = efab.array();
        amrex::Array4<amrex::Real> const& B_loc = bfab.array(); 

        Theta_x(P,particles,n_particles,geom,E_loc,B_loc,dt);
    }
    P.Redistribute();
    P.updateNeighbors();
    E.FillBoundary();
    B.FillBoundary();
}
void G_Theta_y(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab&E, amrex::MultiFab&B,double dt ){
    for (CParIter pti(P, 0); pti.isValid(); ++pti) {
        auto&  particles = pti.GetArrayOfStructs();
        const auto& n_particles = P.GetNeighbors(0,pti.index(),pti.LocalTileIndex());
        amrex::FArrayBox& efab = E[pti];
        amrex::FArrayBox& bfab =B[pti];
        amrex::Array4<amrex::Real> const& E_loc = efab.array();
        amrex::Array4<amrex::Real> const& B_loc = bfab.array(); 

        Theta_y(P,particles,n_particles,geom,E_loc,B_loc,dt);
    }
    P.Redistribute();
    P.updateNeighbors();
    E.FillBoundary();
    B.FillBoundary();
}
void G_Theta_z(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab&E, amrex::MultiFab&B,double dt ){
    for (CParIter pti(P, 0); pti.isValid(); ++pti) {
        auto&  particles = pti.GetArrayOfStructs();
        const auto& n_particles = P.GetNeighbors(0,pti.index(),pti.LocalTileIndex());
        amrex::FArrayBox& efab = E[pti];
        amrex::FArrayBox& bfab =B[pti];
        amrex::Array4<amrex::Real> const& E_loc = efab.array();
        amrex::Array4<amrex::Real> const& B_loc = bfab.array(); 

        Theta_z(P,particles,n_particles,geom,E_loc,B_loc,dt);
    }
    P.Redistribute();
    P.updateNeighbors();
    E.FillBoundary();
    B.FillBoundary();
}
