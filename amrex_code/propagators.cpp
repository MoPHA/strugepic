#include "AMReX_Array.H"
#include "AMReX_Geometry.H"
#include "AMReX_Loop.H"
#include "AMReX_MultiFab.H"
#include "particle_defs.hpp"
#include "propagators.hpp"
#include "amrex_util.hpp"


// Soft sine plane wave source in x-direction
E_source::E_source(amrex::Geometry geom, amrex::MultiFab &E,int pos,int comp,double E0 ,double omega ,double dt) : 
     geom(geom),
     E(E),
     dt(dt),
     omega(omega),
     E0(E0),
     comp(comp),
     pos(pos)
     {}
// Why sin and factor two ?? and not derivative
void E_source::operator()(double t){ 
    for (amrex::MFIter mfi(E); mfi.isValid(); ++mfi){
        amrex::Array4<amrex::Real> const& b = E.array(mfi); 
        const auto box= mfi.validbox();
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i,int j,int k ){
            if(i == pos ){
                b(i,j,k,comp) +=2*E0*sin(omega*t)*dt;
            }
        });
    }
    E.FillBoundary(geom.periodicity());
 }


// Compact hard source sin(wt)
// Keept here for debugging purpouse
void E_source_hard(amrex::Geometry geom,amrex::MultiFab &E,int pos,double omega,double t){
    for (amrex::MFIter mfi(E); mfi.isValid(); ++mfi){
        amrex::Array4<amrex::Real> const& b = E.array(mfi); 
        const auto box= mfi.validbox();
   amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i,int j,int k ){
           if(i == pos ){
            b(i,k,j,Y) =sin(omega*t);
            }
         });
    }
    E.FillBoundary(geom.periodicity());
}


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
    push_ff<-1,1>(geom,bx,E,B,curl_fdiff_1,dt);
}

void push_E_B(const amrex::Geometry geom, amrex::Box const& bx,  amrex::Array4<amrex::Real> const& E, amrex::Array4<amrex::Real const> const& B,double dt){
    push_ff<1,1>(geom,bx,B,E,curl_bdiff_1,dt);
    }/*
void push_B_E(const amrex::Geometry geom, amrex::Box const& bx,  amrex::Array4<amrex::Real> const& B, amrex::Array4<amrex::Real const> const& E,double dt){
   const auto ics = geom.InvCellSize() ;
   const auto domain =geom.Domain();
   const auto lo = amrex::lbound(domain);
   const auto hi = amrex::ubound(domain);
   if (!geom.isPeriodic(X)){
   amrex::ParallelFor(bx,  [=] AMREX_GPU_DEVICE (int i,int j,int k ){
         if(!( (i == lo.x || i==hi.x ) ) ){
         double B_old_x =0;
         double B_old_y =0;
         double B_old_z =0;
         if(i == lo.x+1 || i==hi.x-1){
           B_old_x=B(i,j,k,X); 
           B_old_y=B(i,j,k,Y); 
           B_old_z=B(i,j,k,Z); 
         }
         auto curl = curl_fdiff_1(E,i,j,k,ics);
         B(i,j,k,X) -=dt*curl[0];  
         B(i,j,k,Y) -=dt*curl[1];
         B(i,j,k,Z) -=dt*curl[2];

         if(i== lo.x+1){
         B(i-1,j,k,X) = (1-dt)*B(i-1,j,k,X)+B_old_x*dt;
         B(i-1,j,k,Y) = (1-dt)*B(i-1,j,k,Y)+B_old_y*dt;
         B(i-1,j,k,Z) = (1-dt)*B(i-1,j,k,Z)+B_old_z*dt;
         }
         else if( i== hi.x-1){
         B(i+1,j,k,X) = (1-dt)*B(i+1,j,k,X)+B_old_x*dt;
         B(i+1,j,k,Y) = (1-dt)*B(i+1,j,k,Y)+B_old_y*dt;
         B(i+1,j,k,Z) = (1-dt)*B(i+1,j,k,Z)+B_old_z*dt;
         }
         }
   });
   }
   else{
   amrex::ParallelFor(bx,  [=] AMREX_GPU_DEVICE (int i,int j,int k ){
         auto curl = curl_fdiff_1(E,i,j,k,ics);
         B(i,j,k,X) -=dt*curl[0];  
         B(i,j,k,Y) -=dt*curl[1];
         B(i,j,k,Z) -=dt*curl[2];
   });
   }
}

void push_E_B(const amrex::Geometry geom, amrex::Box const& bx,  amrex::Array4<amrex::Real> const& E, amrex::Array4<amrex::Real const> const& B,double dt){
   const auto ics = geom.InvCellSize() ;
   const auto domain =geom.Domain();
   const auto lo = amrex::lbound(domain);
   const auto hi = amrex::ubound(domain);
   if (!geom.isPeriodic(X)){
   amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i,int j,int k ){
         if(!( (i == lo.x || i==hi.x ) && !geom.isPeriodic(X)) ){
         auto curl = curl_bdiff_1(B,i,j,k,ics);
         double E_old_x=0;
         double E_old_y=0;
         double E_old_z=0;
         if( i== lo.x+1 || i== hi.x -1){
            E_old_x = E(i,j,k,X);
            E_old_y = E(i,j,k,Y);
            E_old_z = E(i,j,k,Z);
         }
         E(i,j,k,X) +=dt*curl[0];  
         E(i,j,k,Y) +=dt*curl[1];
         E(i,j,k,Z) +=dt*curl[2];
         if(i== lo.x+1){
         E(i-1,j,k,X) = (1-dt)*E(i-1,j,k,X)+E_old_x*dt;
         E(i-1,j,k,Y) = (1-dt)*E(i-1,j,k,Y)+E_old_y*dt;
         E(i-1,j,k,Z) = (1-dt)*E(i-1,j,k,Z)+E_old_z*dt;
         }
         else if( i== hi.x-1){
         E(i+1,j,k,X) = (1-dt)*E(i+1,j,k,X)+E_old_x*dt;
         E(i+1,j,k,Y) = (1-dt)*E(i+1,j,k,Y)+E_old_y*dt;
         E(i+1,j,k,Z) = (1-dt)*E(i+1,j,k,Z)+E_old_z*dt;
         }
        } 
         });
   }
   else{
   amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i,int j,int k ){
         auto curl = curl_bdiff_1(B,i,j,k,ics);
         E(i,j,k,X) +=dt*curl[0];  
         E(i,j,k,Y) +=dt*curl[1];
         E(i,j,k,Z) +=dt*curl[2];
         });

   }
}
*/

void Theta_B(const amrex::Geometry geom, amrex::Box const& bx,amrex::Array4<amrex::Real> const& E,amrex::Array4<amrex::Real const> const& B,double dt ){
    push_E_B(geom,bx,E,B,dt);
}




void push_V_E( CParticles&particles, const amrex::Geometry geom,amrex::Array4<amrex::Real const> const& E ,double dt){

    // Basis functions are supported over two cells in each direction
    const int idx_list[4]={-1,0,1,2};
    const auto low =geom.ProbLo();
    const auto Ics = geom.InvCellSize();
    
    for(auto& p : particles){
    const double m= p.rdata(M);
    const double q= p.rdata(Q);
    const double coef = dt*q/m;
        auto coord =get_point_cell(geom,{p.pos(X),p.pos(Y),p.pos(Z)}) ;
       
        double dvx=0;
        double dvy=0;
        double dvz=0;


        std::array<double,4> W1X={0,0,0,0};
        std::array<double,4> W1Y={0,0,0,0};
        std::array<double,4> W1Z={0,0,0,0};
        std::array<double,4> W12X={0,0,0,0};
        std::array<double,4> W12Y={0,0,0,0};
        std::array<double,4> W12Z={0,0,0,0};

        auto nz = (p.pos(Z)-low[Z])*Ics[Z]; 
        for(auto k:idx_list){
            auto cz = coord[Z]+k;
            W1Z[k+1]=W1(nz-cz);
            W12Z[k+1]=Wp(nz-cz);
        }
        auto ny = (p.pos(Y)-low[Y])*Ics[Y]; 
        for(auto j:idx_list){
            auto cy = coord[Y]+j;
            W1Y[j+1]=W1(ny-cy);
            W12Y[j+1]=Wp(ny-cy);
        }
        auto nx = (p.pos(X)-low[X])*Ics[X]; 
        for(auto i:idx_list){
            auto cx = coord[X]+i;
            W1X[i+1]=W1(nx-cx);
            W12X[i+1]=Wp(nx-cx);
        }
        


        for(auto k: idx_list){
            auto cz = coord[Z]+k;
            for(auto j: idx_list){
                  auto cy = coord[Y]+j;
                for(auto i: idx_list){
                  auto cx = coord[X]+i;
                   dvx+= E(cx,cy,cz,X)*W12X[i+1]*W1Y[j+1]*W1Z[k+1]; 
                   dvy+= E(cx,cy,cz,Y)*W1X[i+1]*W12Y[j+1]*W1Z[k+1];
                   dvz+= E(cx,cy,cz,Z)*W1X[i+1]*W1Y[j+1]*W12Z[k+1];
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

    B.FillBoundary(geom.periodicity());


}



 void G_Theta_B(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab &E, amrex::MultiFab &B,double dt ){
    
    for (amrex::MFIter mfi(E); mfi.isValid(); ++mfi){
        const amrex::Box& box = mfi.validbox();
        auto const& E_loc = E.array(mfi);
        auto const& B_loc = B.const_array(mfi); 
        Theta_B(geom,box,E_loc,B_loc,dt);

    }

    E.FillBoundary(geom.periodicity());

}
