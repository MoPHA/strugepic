#include "AMReX_Array.H"
#include "AMReX_Geometry.H"
#include "AMReX_Loop.H"
#include "AMReX_MultiFab.H"
#include "strugepic_defs.hpp"
#include "strugepic_util.hpp"
#include "strugepic_propagators.hpp"




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
void E_source::operator()(double t){ 
    for (amrex::MFIter mfi(E); mfi.isValid(); ++mfi){
        amrex::Array4<amrex::Real> const& b = E.array(mfi); 
        const auto box= mfi.validbox();
        amrex::ParallelFor(box, [=]  (int i,int j,int k ){
            if(i == pos ){
                b(i,j,k,comp) +=2*E0*sin(omega*t)*dt;
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
inline std::array<amrex::Real,3> curl_fdiff_1(amrex::Array4<amrex::Real const> const& a ,int i , int j ,int k, const double* ics){
    return {
        ((a(i,j+1,k,Z)-a(i,j,k,Z))*ics[Y]-( a(i,j,k+1,Y)-a(i,j,k,Y))*ics[Z]),
        ((a(i,j,k+1,X)-a(i,j,k,X))*ics[Z]-( a(i+1,j,k,Z)-a(i,j,k,Z))*ics[X]),
        ((a(i+1,j,k,Y)-a(i,j,k,Y))*ics[X]-( a(i,j+1,k,X)-a(i,j,k,X))*ics[Y])
    };
}

// backward difference, first order
inline std::array<amrex::Real,3> curl_bdiff_1(amrex::Array4<amrex::Real const> const& a ,int i , int j ,int k, const double* ics){
    return {
        ((a(i,j,k,Z)-a(i,j-1,k,Z))*ics[Y]-( a(i,j,k,Y)-a(i,j,k-1,Y))*ics[Z]),
        ((a(i,j,k,X)-a(i,j,k-1,X))*ics[Z]-( a(i,j,k,Z)-a(i-1,j,k,Z))*ics[X]),
        ((a(i,j,k,Y)-a(i-1,j,k,Y))*ics[X]-( a(i,j,k,X)-a(i,j-1,k,X))*ics[Y])
    };
}

void E_curl(const amrex::Geometry geom, amrex::Box const& bx,  amrex::Array4<amrex::Real const> const& E, amrex::Array4<amrex::Real > const& B,double dt){

   amrex::ParallelFor(bx, [=]  (int i,int j,int k ){
         auto res = curl_fdiff_1(E,i,j,k,geom.InvCellSize());
         B(i,j,k,X) =B(i,j,k,X)-dt*res[0] ;  
         B(i,j,k,Y) =B(i,j,k,Y)-dt*res[1] ;
         B(i,j,k,Z) =B(i,j,k,Z)-dt*res[2] ;   
        }); 
}

void B_curl(const amrex::Geometry geom, amrex::Box const& bx,  amrex::Array4<amrex::Real const> const& B, amrex::Array4<amrex::Real> const& E,double dt){
   amrex::ParallelFor(bx, [=]  (int i,int j,int k ){
         auto res = curl_bdiff_1(B,i,j,k,geom.InvCellSize());
         E(i,j,k,X) =E(i,j,k,X)+dt*res[0] ;  
         E(i,j,k,Y) =E(i,j,k,Y)+dt*res[1] ;
         E(i,j,k,Z) =E(i,j,k,Z)+dt*res[2] ;   
        }); 

}

void push_B_E(const amrex::Geometry geom, amrex::Box const& bx,  amrex::Array4<amrex::Real> const& B, amrex::Array4<amrex::Real const> const& E,double dt){
    push_ff<E_curl,MABC<X>>(geom,bx,E,B,dt);
}

void push_E_B(const amrex::Geometry geom, amrex::Box const& bx,  amrex::Array4<amrex::Real> const& E, amrex::Array4<amrex::Real const> const& B,double dt){
    push_ff<B_curl,MABC<X>>(geom,bx,B,E,dt);
    }


 void G_Theta_B(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab &E, amrex::MultiFab &B,double dt ){
    
    B.FillBoundary(geom.periodicity());
    for (amrex::MFIter mfi(E); mfi.isValid(); ++mfi){
        const amrex::Box& box = mfi.validbox();
        auto const& E_loc = E.array(mfi);
        auto const& B_loc = B.const_array(mfi);  
        push_E_B(geom,box,E_loc,B_loc,dt);
    }


}
