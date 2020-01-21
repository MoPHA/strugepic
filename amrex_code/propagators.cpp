#include "AMReX_Array.H"
#include "AMReX_Geometry.H"
#include "AMReX_Loop.H"
#include "particle_defs.hpp"
#include "propagators.hpp"
#include "amrex_util.hpp"


// Testing another thing 



// Central difference first order
 std::array<amrex::Real,3> curl_cdiff_1(amrex::Array4<amrex::Real const> const& a ,int i , int j ,int k, const double* ics){
    return {
        ((a(i,j+1,k,Z)-a(i,j-1,k,Z))*ics[Y]*0.5-( a(i,j,k+1,Y)-a(i,j,k-1,Y))*ics[Z]*0.5),
        ((a(i,j,k+1,X)-a(i,j,k-1,X))*ics[Z]*0.5-( a(i+1,j,k,Z)-a(i-1,j,k,Z))*ics[X]*0.5),    
        ((a(i+1,j,k,Y)-a(i-1,j,k,Y))*ics[X]*0.5-( a(i,j+1,k,X)-a(i,j-1,k,X))*ics[Y]*0.5)
    };
}


// forward difference, first order
 std::array<amrex::Real,3> curl_fdiff_1(amrex::Array4<amrex::Real const> const& a ,int i , int j ,int k, const double* ics){
    return {
        ((a(i,j+1,k,Z)-a(i,j,k,Z))*ics[Y]-( a(i,j,k+1,Y)-a(i,j,k,Y))*ics[Z]),
        ((a(i,j,k+1,X)-a(i,j,k,X))*ics[Z]-( a(i+1,j,k,Z)-a(i,j,k,Z))*ics[X]),
        ((a(i+1,j,k,Y)-a(i,j,k,Y))*ics[X]-( a(i,j+1,k,X)-a(i,j,k,X))*ics[Y])
    };
}

// backward difference, first order
 std::array<amrex::Real,3> curl_bdiff_1(amrex::Array4<amrex::Real const> const& a ,int i , int j ,int k, const double* ics){
    return {
        ((a(i,j,k,Z)-a(i,j-1,k,Z))*ics[Y]-( a(i,j,k,Y)-a(i,j,k-1,Y))*ics[Z]),
        ((a(i,j,k,X)-a(i,j,k-1,X))*ics[Z]-( a(i,j,k,Z)-a(i-1,j,k,Z))*ics[X]),
        ((a(i,j,k,Y)-a(i-1,j,k,Y))*ics[X]-( a(i,j,k,X)-a(i,j-1,k,X))*ics[Y])
    };
}


void push_B_E(const amrex::Geometry geom, amrex::Box const& bx,  amrex::Array4<amrex::Real> const& B, amrex::Array4<amrex::Real const> const& E,double dt){
   const auto ics = geom.InvCellSize() ;
   amrex::ParallelFor(bx,  [=] AMREX_GPU_DEVICE (int i,int j,int k ){
         auto curl = curl_fdiff_1(E,i,j,k,ics);
         B(i,j,k,X) -=dt*curl[0];  
         B(i,j,k,Y) -=dt*curl[1];
         B(i,j,k,Z) -=dt*curl[2];
         });

}

void push_E_B(const amrex::Geometry geom, amrex::Box const& bx,  amrex::Array4<amrex::Real> const& E, amrex::Array4<amrex::Real const> const& B,double dt){
   const auto ics = geom.InvCellSize() ;
   const auto lo = amrex::lbound(bx);
   const auto hi = amrex::ubound(bx);
   amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i,int j,int k ){
         if( (i == lo.x || i==hi.x ) && !geom.isPeriodic(X) ){
         return;
         }
        auto curl = curl_bdiff_1(B,i,j,k,ics);
         E(i,j,k,X) +=dt*curl[0];  
         E(i,j,k,Y) +=dt*curl[1];
         E(i,j,k,Z) +=dt*curl[2];
         });
}


void Theta_B(const amrex::Geometry geom, amrex::Box const& bx,amrex::Array4<amrex::Real> const& E,amrex::Array4<amrex::Real const> const& B,double dt ){
    push_E_B(geom,bx,E,B,dt);
}




void push_V_E( CParticles&particles, const amrex::Geometry geom,amrex::Array4<amrex::Real const> const& E ,double dt){

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
        amrex::Array4<amrex::Real const> const& E_loc = E.const_array(pti);
         
        push_V_E(particles,geom,E_loc,dt);
    }

    for (amrex::MFIter mfi(E); mfi.isValid(); ++mfi){
        const amrex::Box& box = mfi.validbox();
        amrex::Array4<amrex::Real const> const& E_loc = E.const_array(mfi);
        amrex::Array4<amrex::Real> const& B_loc = B.array(mfi); 
        push_B_E(geom,box, B_loc,E_loc,dt);

    }

    P.updateNeighbors();
    B.FillBoundary(geom.periodicity());


}



 void G_Theta_B(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab &E, amrex::MultiFab &B,double dt ){
    
    for (amrex::MFIter mfi(E); mfi.isValid(); ++mfi){
        const amrex::Box& box = mfi.validbox();
        auto const& E_loc = E.array(mfi);
        auto const& B_loc = B.const_array(mfi); 
        Theta_B(geom,box,E_loc,B_loc,dt);

    }

    P.updateNeighbors();
    E.FillBoundary(geom.periodicity());

}
