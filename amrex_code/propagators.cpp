#include "AMReX_Array.H"
#include "AMReX_Geometry.H"
#include "particle_defs.hpp"
#include "propagators.hpp"
#include "amrex_util.hpp"



void Theta_B(amrex::Box const& bx,amrex::Array4<amrex::Real> const& E,amrex::Array4<amrex::Real> const& B,double dt ){
    push_E_B(bx,E,B,dt);
}


void push_B_E( amrex::Box const& bx,  amrex::Array4<amrex::Real> const& B, amrex::Array4<amrex::Real> const& E,double dt){
   const auto lo = amrex::lbound(bx);
   const auto hi = amrex::ubound(bx);
   for     (int k = lo.z; k <= hi.z; ++k) {
     for   (int j = lo.y; j <= hi.y; ++j) {
       for (int i = lo.x; i <= hi.x; ++i) { 
         B(i,j,k,X) -=dt*((E(i,j+1,k,Z)-E(i,j,k,Z))-( E(i,j,k+1,Y)-E(i,j,k,Y)));  
         B(i,j,k,Y) -=dt*((E(i,j,k+1,X)-E(i,j,k,X))-( E(i+1,j,k,Z)-E(i,j,k,Z)));
         B(i,j,k,Z) -=dt*((E(i+1,j,k,Y)-E(i,j,k,Y))-( E(i,j+1,k,X)-E(i,j,k,X)));
       }
     }
   }

}

void push_E_B( amrex::Box const& bx,  amrex::Array4<amrex::Real> const& E, amrex::Array4<amrex::Real> const& B,double dt){
   const auto lo = amrex::lbound(bx);
   const auto hi = amrex::ubound(bx);
   for     (int k = lo.z; k <= hi.z; ++k) {
     for   (int j = lo.y; j <= hi.y; ++j) {
       for (int i = lo.x; i <= hi.x; ++i) { 
         E(i,j,k,X) +=dt*((B(i,j+1,k,Z)-B(i,j,k,Z))-( B(i,j,k+1,Y)-B(i,j,k,Y)));  
         E(i,j,k,Y) +=dt*((B(i,j,k+1,X)-B(i,j,k,X))-( B(i+1,j,k,Z)-B(i,j,k,Z)));
         E(i,j,k,Z) +=dt*((B(i+1,j,k,Y)-B(i,j,k,Y))-( B(i,j+1,k,X)-B(i,j,k,X)));
       }
     }
   }

}


void push_V_E( CParticles&particles, const amrex::Geometry geom,amrex::Array4<amrex::Real> const& E ,double dt){

    // Basis functions are supported over two cells in each direction
    const int idx_list[4]={-1,0,1,2};
    // All the particles in a give container have the same mass and charge
    // We could probably do without saving the mass and charge for each particle
    const double m= particles[0].rdata(M);
    const double q= particles[0].rdata(Q);
    const double coef = dt*q/m;
    const auto low =geom.ProbLo();
    const auto Ics = geom.InvCellSize();
    
    for(auto& p : particles){
        auto coord =get_point_cell(geom,{p.pos(X),p.pos(Y),p.pos(Z)}) ;
       
        double dvx=0;
        double dvy=0;
        double dvz=0;

        for(auto k: idx_list){
            for(auto j: idx_list){
                for(auto i: idx_list){
                  auto cx = coord[X]+i;
                  auto cy = coord[Y]+j;
                  auto cz = coord[Z]+k;
                  using namespace strugepic;
                   dvx+= E(cx,cy,cz,X)*W12((p.pos(X)-low[X])*Ics[X]-cx)*W1((p.pos(Y)-low[Y])*Ics[Y]-cy)*W1((p.pos(Z)-low[Z])*Ics[Z]-cz); 
                   dvy+= E(cx,cy,cz,Y)*W1((p.pos(X)-low[X])*Ics[X]-cx)*W12((p.pos(Y)-low[Y])*Ics[Y]-cy)*W1((p.pos(Z)-low[Z])*Ics[Z]-cz);
                   dvz+= E(cx,cy,cz,Z)*W1((p.pos(X)-low[X])*Ics[X]-cx)*W1((p.pos(Y)-low[Y])*Ics[Y]-cy)*W12((p.pos(Z)-low[Z])*Ics[Z]-cz);
                }
            }
        }
        p.rdata(VX)+=dvx*coef;
        p.rdata(VY)+=dvy*coef;
        p.rdata(VZ)+=dvz*coef;
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
        push_B_E(box, B_loc,E_loc,dt);

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
    E.FillBoundary(geom.periodicity());

}
