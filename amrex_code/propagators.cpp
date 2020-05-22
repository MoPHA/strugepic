#include "AMReX_Array.H"
#include "AMReX_Geometry.H"
#include "AMReX_Loop.H"
#include "AMReX_MultiFab.H"
#include "particle_defs.hpp"
#include "propagators.hpp"
#include "amrex_util.hpp"


void fill_extra_halo(amrex::Geometry geom, amrex::Array4<amrex::Real> const& A,amrex::Box bx,int ng){
    for(auto coord: {X,Y,Z}){
        auto coord_u = (coord+1)%3;
        auto coord_l = (coord+2)%3;
        const auto domainLo = geom.Domain().loVect();
        const auto subgridLo =bx.loVect();
        const auto domainHi = geom.Domain().hiVect();
        const auto subgridHi = bx.hiVect();
        auto missing_L =  ng-(domainHi[coord]-domainLo[coord]+1);
        if(domainLo[coord] == subgridLo[coord] && domainHi[coord]==subgridHi[coord] && geom.isPeriodic(coord)
                && missing_L > 0){
                const auto e_low=bx.loVect();
                const auto e_high=bx.hiVect();
                auto elx=e_high[X]-e_low[X] +1;
                auto ely=e_high[Y]-e_low[Y] +1;
                auto elz=e_high[Z]-e_low[Z] +1;
        
             for(int u=0; u<=(subgridHi[coord_u]-subgridLo[coord_u]);u++){
                for(int l=0; l<=(subgridHi[coord_l]-subgridLo[coord_l]);l++){
                   for(int c =0; c < missing_L ; c++){
                        
                        std::array<int,3> ext_coord;
                        ext_coord[coord_u] = u+subgridLo[coord_u];
                        ext_coord[coord_l] = l+subgridLo[coord_l];
                        ext_coord[coord] = c+subgridHi[coord]+1+(subgridHi[coord]-subgridLo[coord]+1);
                    std::array<int,3> int_coord;
                    int_coord[coord_l]=ext_coord[coord_l];
                    int_coord[coord_u]=ext_coord[coord_u];
                    int_coord[coord]=ext_coord[coord];
                    int_coord[X]=mod((int_coord[X]-e_low[X]), elx)+e_low[X];
                    int_coord[Y]=mod((int_coord[Y]-e_low[Y]), ely)+e_low[Y];
                    int_coord[Z]=mod((int_coord[Z]-e_low[Z]), elz)+e_low[Z];
                A(ext_coord[X],ext_coord[Y],ext_coord[Z],X)=A(int_coord[X],int_coord[Y],int_coord[Z],X);
                A(ext_coord[X],ext_coord[Y],ext_coord[Z],Y)=A(int_coord[X],int_coord[Y],int_coord[Z],Y);
                A(ext_coord[X],ext_coord[Y],ext_coord[Z],Z)=A(int_coord[X],int_coord[Y],int_coord[Z],Z);
                   }
                }
             }
             for(int u=0; u<=(subgridHi[coord_u]-subgridLo[coord_u]);u++){
                for(int l=0; l<=(subgridHi[coord_l]-subgridLo[coord_l]);l++){
                   for(int c =0; c < missing_L ; c++){
                        
                        std::array<int,3> ext_coord;
                        ext_coord[coord_u] = u+subgridLo[coord_u];
                        ext_coord[coord_l] = l+subgridLo[coord_l];
                        ext_coord[coord] = c+subgridLo[coord]-ng;
                    std::array<int,3> int_coord;
                    int_coord[coord_l]=ext_coord[coord_l];
                    int_coord[coord_u]=ext_coord[coord_u];
                    int_coord[coord]=ext_coord[coord];
                    int_coord[X]=mod((int_coord[X]-e_low[X]), elx)+e_low[X];
                    int_coord[Y]=mod((int_coord[Y]-e_low[Y]), ely)+e_low[Y];
                    int_coord[Z]=mod((int_coord[Z]-e_low[Z]), elz)+e_low[Z];
                A(ext_coord[X],ext_coord[Y],ext_coord[Z],X)=A(int_coord[X],int_coord[Y],int_coord[Z],X);
                A(ext_coord[X],ext_coord[Y],ext_coord[Z],Y)=A(int_coord[X],int_coord[Y],int_coord[Z],Y);
                A(ext_coord[X],ext_coord[Y],ext_coord[Z],Z)=A(int_coord[X],int_coord[Y],int_coord[Z],Z);
                   }
                }
             }

        }
    }    
}


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
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i,int j,int k ){
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

   amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i,int j,int k ){
         auto res = curl_fdiff_1(E,i,j,k,geom.InvCellSize());
         B(i,j,k,X) =B(i,j,k,X)-dt*res[0] ;  
         B(i,j,k,Y) =B(i,j,k,Y)-dt*res[1] ;
         B(i,j,k,Z) =B(i,j,k,Z)-dt*res[2] ;   
        }); 
}

void B_curl(const amrex::Geometry geom, amrex::Box const& bx,  amrex::Array4<amrex::Real const> const& B, amrex::Array4<amrex::Real> const& E,double dt){
   amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i,int j,int k ){
         auto res = curl_bdiff_1(B,i,j,k,geom.InvCellSize());
         E(i,j,k,X) =E(i,j,k,X)+dt*res[0] ;  
         E(i,j,k,Y) =E(i,j,k,Y)+dt*res[1] ;
         E(i,j,k,Z) =E(i,j,k,Z)+dt*res[2] ;   
        }); 

}

void push_B_E(const amrex::Geometry geom, amrex::Box const& bx,  amrex::Array4<amrex::Real> const& B, amrex::Array4<amrex::Real const> const& E,double dt){
    push_ff<1>(geom,bx,E,B,E_curl,MABC<X>,dt);
}

void push_E_B(const amrex::Geometry geom, amrex::Box const& bx,  amrex::Array4<amrex::Real> const& E, amrex::Array4<amrex::Real const> const& B,double dt){
    push_ff<1>(geom,bx,B,E,B_curl,MABC<X>,dt);
    }



void G_Theta_E(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab&E, amrex::MultiFab&B,double dt ){


    for (CParIter pti(P, 0); pti.isValid(); ++pti) {
        auto&  particles = pti.GetArrayOfStructs();
        amrex::Array4<amrex::Real const> const& E_loc = E.const_array(pti);
         
        push_V_E<WRANGE>(particles,geom,E_loc,dt);
    }

    for (amrex::MFIter mfi(E); mfi.isValid(); ++mfi){
        const amrex::Box& box = mfi.validbox();
        amrex::Array4<amrex::Real const> const& E_loc = E.const_array(mfi);
        amrex::Array4<amrex::Real> const& B_loc = B.array(mfi); 
        push_B_E(geom,box, B_loc,E_loc,dt);
        fill_extra_halo(geom,B_loc,box,B.nGrow());
    }

    B.FillBoundary(geom.periodicity());

}


 void G_Theta_B(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab &E, amrex::MultiFab &B,double dt ){
    
    for (amrex::MFIter mfi(E); mfi.isValid(); ++mfi){
        const amrex::Box& box = mfi.validbox();
        auto const& E_loc = E.array(mfi);
        auto const& B_loc = B.const_array(mfi);  
        push_E_B(geom,box,E_loc,B_loc,dt);
        fill_extra_halo(geom,E_loc,box,B.nGrow());
    }

    E.FillBoundary(geom.periodicity());

}
