#ifndef PROPAGATORS
#define PROPAGATORS
#include<AMReX_Geometry.H>
#include <array>
#include "strugepic_w.hpp"
#include "AMReX_Box.H"
#include "AMReX_BoxArray.H"
#include "AMReX_BoxDomain.H"
#include "AMReX_MultiFab.H"
#include "strugepic_defs.hpp"
#include "strugepic_util.hpp"
#include<cmath>






class E_source
{
    public:
        E_source(amrex::Geometry , amrex::MultiFab &E,int pos,int comp,double E0,double omega ,double dt); 
        void operator()(double t);
    private:
        const amrex::Geometry geom;
        amrex::MultiFab &E;
        const int comp;
        const double dt;
        const double omega;
        const int pos;
        const double E0;
};


// These are All local update functions, I.e they operate only on local data 

void push_B_E(const amrex::Geometry geom, amrex::Box const& bx,  amrex::Array4<amrex::Real> const& B, amrex::Array4<amrex::Real const> const& E,double dt);

void push_E_B(const amrex::Geometry geom, amrex::Box const& bx,  amrex::Array4<amrex::Real> const& E, amrex::Array4<amrex::Real const> const& B,double dt);

void Theta_E(const amrex::Geometry geom,amrex::Box const& bx,amrex::Array4<amrex::Real const> const& E,amrex::Array4<amrex::Real> const& B,CParticles&particles,double dt );

void Theta_B(const amrex::Geometry geom,amrex::Box const& bx,amrex::Array4<amrex::Real> const& E,amrex::Array4<amrex::Real const> const& B,double dt );

template<int W_range>
void push_V_E( CParticles&particles, const amrex::Geometry geom,amrex::Array4<amrex::Real const> const& E ,double dt);

// Global update, they also handle triggering global communication
void G_Theta_B(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab &E, amrex::MultiFab &B,double dt );


template<int W_range>
void G_Theta_E(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab&E, amrex::MultiFab&B,double dt ){


    E.FillBoundary(geom.periodicity());
    for (CParIter pti(P, 0); pti.isValid(); ++pti) {
        auto&  particles = pti.GetArrayOfStructs();
        amrex::Array4<amrex::Real const> const& E_loc = E.const_array(pti);
        push_V_E<W_range>(particles,geom,E_loc,dt);
    }
    for (amrex::MFIter mfi(E); mfi.isValid(); ++mfi){
        const amrex::Box& box = mfi.validbox();
        amrex::Array4<amrex::Real const> const& E_loc = E.const_array(mfi);
        amrex::Array4<amrex::Real> const& B_loc = B.array(mfi); 
        push_B_E(geom,box, B_loc,E_loc,dt);
    }


}


typedef bool (*SEG_BOUNDARY)(const amrex::Geometry,amrex::Real *, int *);
typedef void (*PART_BOUNDARY)(CParticle const&, amrex::Real * );

// This is for one particle type
// If there are several you need to do this again
// 0 -> x  , 1-> y 2->z
template<int comp,int W_range,SEG_BOUNDARY F_SEG ,PART_BOUNDARY F_PART>
void Theta(CParticles&particles, const amrex::Geometry geom,amrex::Array4<amrex::Real> const& E ,amrex::Array4<amrex::Real> const& B  ,amrex::Box bx,double dt){ 
    const auto _Ics = geom.InvCellSize();
    const long np = particles.numParticles();
    


    constexpr int W1_r=W_range*2;
    constexpr int W1_li=-W_range+1;
    constexpr int W1_hi= W_range;
    constexpr int Wp_li= W1_li;
    constexpr int Wp_hi= W1_hi-1;
    constexpr int Wp_r=W1_r-1;

    auto comp_u = (comp+1)%3;
    auto comp_l = (comp+2)%3;
    auto comp_Cs = geom.CellSize(comp);
    amrex::GpuArray<float,3> lb;
    amrex::GpuArray<float,3> Ics;
    lb[X] =geom.ProbLo(X);
    lb[Y] =geom.ProbLo(Y);
    lb[Z] =geom.ProbLo(Z);
    Ics[X]=_Ics[X];
    Ics[Y]=_Ics[Y];
    Ics[Z]=_Ics[Z];

    // Remote particle grid lower corner
    amrex::ParallelFor(np,          
            [=] AMREX_GPU_DEVICE (long i)
            {
                    auto &p = particles[i];

    const double m= p.rdata(M);
    const double q= p.rdata(Q);
    const  double B_coef = q/m*comp_Cs;
    const  double E_coef = q*Ics[X]*Ics[Y]*Ics[Z]*comp_Cs;
        double new_pos=p.pos(comp)+dt*p.rdata(comp+2);
        int coord[3];
        coord[X]=floor((p.pos(X) -lb[X])*Ics[X]);
        coord[Y]=floor((p.pos(Y) -lb[Y])*Ics[Y]);
        coord[Z]=floor((p.pos(Z) -lb[Z])*Ics[Z]);
        

        amrex::Real res_c1=0;
        amrex::Real res_c2=0;


        
            amrex::Array1D<amrex::Real,0,W1_r> comp_uW1={0};
            amrex::Array1D<amrex::Real,0,W1_r> comp_lW1={0};
            amrex::Array1D<amrex::Real,0,W1_r> comp_uWp={0};
            amrex::Array1D<amrex::Real,0,W1_r> comp_lWp={0};




            auto nl=(p.pos(comp_l)-lb[comp_l])*Ics[comp_l];
            int idx=0;
            for(int l=W1_li; l<=W1_hi;l++){
                auto cl=coord[comp_l]+(l); 
                comp_lW1[idx]=W1(nl-cl);
                idx++;

            }idx=0;
            for(int l=Wp_li; l<=Wp_hi;l++){
                auto cl=coord[comp_l]+(l); 
                comp_lWp[idx]=Wp(nl-cl);
                idx++;

            }idx=0;
            auto nu=(p.pos(comp_u)-lb[comp_u])*Ics[comp_u];
            for(int  u=W1_li; u<=W1_hi;u++){
                auto cu=coord[comp_u]+(u); 
                comp_uW1[idx]=W1(nu-cu);
                idx++;
            }idx=0;
        

            for(int  u=Wp_li; u<=Wp_hi;u++){
                auto cu=coord[comp_u]+(u); 
                comp_uWp[idx]=Wp(nu-cu);
                idx++;
            }
        
         amrex::Real seg_points[3];
         int seg_idx[2]; 
         int num_segments=construct_segments(p.pos(comp),new_pos,seg_points,seg_idx);
         bool out_not_periodic=F_SEG(geom,seg_points,seg_idx);

        for(int seg=0;seg<num_segments;seg++){
            coord[comp] = seg_idx[seg];
            auto i_s = seg_points[seg]; 
            auto i_e = seg_points[seg+1];
            std::array<double,Wp_r> compI_W12={0};
           
           auto ncs=(i_s-lb[comp])*Ics[comp];
           auto nce=(i_e-lb[comp])*Ics[comp];
           idx=0;
            for(int c=Wp_li;c<=Wp_hi;c++){
                auto cc=coord[comp]+(c);
                compI_W12[idx]=I_Wp(ncs- cc ,nce- cc);
                idx++;

            }


            int idl=0;
            for(int l=W1_li;l<=W1_hi; l++ ){
                int idu=0;
                for(int u=W1_li;u<=W1_hi;u++){
                   int  idc=0;
                double mul=comp_lW1[idl]*comp_uW1[idu];
                double mulu12l1=comp_uWp[idu]*comp_lW1[idl];
                double mull12u1=comp_lWp[idl]*comp_uW1[idu];
                for(int c=Wp_li;c<=Wp_hi;c++){
                    int cx,cy,cz;

                    if(comp==0){
                     cx=coord[X]+(c);
                     cy=coord[Y]+(u);
                     cz=coord[Z]+(l);
                    }else if(comp==1){
                     cx=coord[X]+(l);
                     cy=coord[Y]+(c);
                     cz=coord[Z]+(u);

                    }
                    else{
                     cx=coord[X]+(u);
                     cy=coord[Y]+(l);
                     cz=coord[Z]+(c);
                    }
                    amrex::Gpu::Atomic::Add(E(cx,cy,cz,comp),-E_coef*mul*compI_W12[idc]);
                    res_c1+=B(cx,cy,cz,comp_u)*mull12u1*compI_W12[idc];
                    res_c2-=B(cx,cy,cz,comp_l)*compI_W12[idc]*mulu12l1;
                    idc++;
                }
                idu++;
            }
                idl++;
            }


        
        }

        // Update position
            if(out_not_periodic){
                F_PART(p,seg_points);                
                //particle_reflect<comp>(p,seg_points);
            }

        // .Redistribute should wrap the particles if periodic
            else{
            p.pos(comp)+=dt*p.rdata(2+comp);
            }
        // B Vel update
        p.rdata( (comp +2 )% 3 +2   )+=B_coef*res_c1;
        p.rdata( (comp+1) % 3 +2  )+=B_coef*res_c2; 

    });
}


template<int W_range>
void push_V_E( CParticles&particles, const amrex::Geometry geom,amrex::Array4<amrex::Real const> const& E ,double dt){

    // Basis functions are supported over two cells in each direction
    const auto low =geom.ProbLo();
    const auto Ics = geom.InvCellSize();
    
    for(auto& p : particles){
    const double m= p.rdata(M);
    const double q= p.rdata(Q);
    const double coef = dt*q/m;
    int coord[3];
    coord[X]=floor((p.pos(X) -low[X])*Ics[X]);
    coord[Y]=floor((p.pos(Y) -low[Y])*Ics[Y]);
    coord[Z]=floor((p.pos(Z) -low[Z])*Ics[Z]);
       
        double dvx=0;
        double dvy=0;
        double dvz=0;

        constexpr int W1_r=W_range*2;
        constexpr int W1_li=-W_range+1;
        constexpr int W1_hi= W_range;
        constexpr int Wp_li= W1_li;
        constexpr int Wp_hi= W1_hi-1;

        std::array<double,W1_r> W1X={0};
        std::array<double,W1_r> W1Y={0};
        std::array<double,W1_r> W1Z={0};
        std::array<double,W1_r> WpX={0};
        std::array<double,W1_r> WpY={0};
        std::array<double,W1_r> WpZ={0};

        // Particle coordinate shifted to cell
        auto nz = (p.pos(Z)-low[Z])*Ics[Z]; 
        auto ny = (p.pos(Y)-low[Y])*Ics[Y]; 
        auto nx = (p.pos(X)-low[X])*Ics[X]; 
      
            int idx=0;
            for(int  i=W1_li; i<=W1_hi;i++){
                auto cx = coord[X]+i;
                W1X[idx]=W1(nx-cx); 
                auto cy = coord[Y]+i;
                W1Y[idx]=W1(ny-cy);
                auto cz = coord[Z]+i;
                W1Z[idx]=W1(nz-cz);
                idx++;
            }
            idx=0;
            for(int i=Wp_li; i<=Wp_hi;i++){
                auto cx = coord[X]+i;
                WpX[idx]=Wp(nx-cx); 
                auto cy = coord[Y]+i;
                WpY[idx]=Wp(ny-cy);
                auto cz = coord[Z]+i;
                WpZ[idx]=Wp(nz-cz);
                idx++;
            
            }

        


            int idk=0;
        for(int  k=W1_li; k<=W1_hi;k++){
            auto cz = coord[Z]+k;
            int idj=0;
            for(int  j=W1_li; j<=W1_hi;j++){
                int idi=0;
                  auto cy = coord[Y]+j;
            for(int  i=W1_li; i<=W1_hi;i++){
                  auto cx = coord[X]+i;
                   dvx+= E(cx,cy,cz,X)*WpX[idi]*W1Y[idj]*W1Z[idk]; 
                   dvy+= E(cx,cy,cz,Y)*W1X[idi]*WpY[idj]*W1Z[idk];
                   dvz+= E(cx,cy,cz,Z)*W1X[idi]*W1Y[idj]*WpZ[idk];
                   idi++;
                }
                idj++;
            }
            idk++;
        }
        p.rdata(VX)+=dvx*coef;
        p.rdata(VY)+=dvy*coef;
        p.rdata(VZ)+=dvz*coef;
        
    }
}



template <int comp,int W_range>
void G_Theta(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab &E, amrex::MultiFab &B,double dt ){

    B.FillBoundary(geom.periodicity());
    E.setBndry(0);
    E.setDomainBndry(0,geom);
    for (amrex::MFIter mfi(E); mfi.isValid(); ++mfi){
        auto box=E.box(mfi.index());
         
        // Each grid,tile has a their own local particle container
        auto& Part = P.GetParticles(0)[std::make_pair(mfi.index(),mfi.LocalTileIndex())];
        auto&  particles = Part.GetArrayOfStructs();
        amrex::FArrayBox& bfab =B[mfi];
        
        amrex::Array4<amrex::Real> const& B_loc = bfab.array(); 
        amrex::Array4<amrex::Real> const& E_loc = E[mfi].array(); 
        Theta<comp,W_range,segment_reflect<comp,W_range>,particle_reflect<comp>>(particles,geom,E_loc,B_loc,box,dt);
    }
    
    E.SumBoundary(geom.periodicity());
    P.Redistribute();
    

    
}


template<int comp,int boundary_size>
void construct_exterior_interior(const amrex::Geometry geom ,std::vector<std::array<int,3>> &exterior,amrex::Box &interior){
   const auto domain=geom.Domain();
   const auto Lo = domain.loVect();
   const auto Hi = domain.hiVect();
   const auto local_lo=interior.loVect();
   const auto local_hi=interior.hiVect();
   if(!geom.isPeriodic(comp) && boundary_size > 0){
        std::array<int,3> point;
        const int up_c =(comp+1) % 3;
        const int lo_c =(comp+2) % 3;
       if(interior.loVect()[comp]==Lo[comp]){
           interior.growLo(comp,-boundary_size);
            for(int k=local_lo[up_c]; k <= local_hi[up_c];k++){
                for(int j=local_lo[lo_c]; j<=local_hi[lo_c];j++ ){
                    for(int i=Lo[comp]; i < Lo[comp]+boundary_size;i++){
                        point[up_c]=k;
                        point[lo_c]=j;
                        point[comp]=i;
                        exterior.push_back(point);

                    }
                }
            }
       }
        
       if(interior.hiVect()[comp]==Hi[comp]){
           interior.growHi(comp,-boundary_size);
            for(int k=local_lo[up_c]; k <= local_hi[up_c];k++){
                for(int j=local_lo[lo_c]; j<=local_hi[lo_c];j++ ){
                    for(int i=Hi[comp]; i > Hi[comp]-boundary_size;i--){
                        point[up_c]=k;
                        point[lo_c]=j;
                        point[comp]=i;
                        exterior.push_back(point);
                    }
                }
            }
       }
   }     
}

template<int comp>
void MABC(const amrex::Geometry geom,amrex::Array4<amrex::Real> const& A,std::vector<std::array<int,3>> const&exterior,double dt){
   const auto domain=geom.Domain();
   const auto Lo = domain.loVect();
   const auto Hi = domain.hiVect();
        for(auto point : exterior){
          if(point[comp] == Lo[comp]){
            const int i = point[X];
            const int j = point[Y];
            const int k = point[Z];
            auto npoint=point;
            npoint[comp]+=1;
           A(i,j,k,X) = (1-dt)*A(i,j,k,X)+A(npoint[X],npoint[Y],npoint[Z],X)*dt;
           A(i,j,k,Y) = (1-dt)*A(i,j,k,Y)+A(npoint[X],npoint[Y],npoint[Z],Y)*dt;
           A(i,j,k,Z) = (1-dt)*A(i,j,k,Z)+A(npoint[X],npoint[Y],npoint[Z],Z)*dt; 
          }
          else if(point[comp] == Hi[comp]){ 
            const int i = point[X];
            const int j = point[Y];
            const int k = point[Z];
            auto npoint=point;
            npoint[comp]-=1;
           A(i,j,k,X) = (1-dt)*A(i,j,k,X)+A(npoint[X],npoint[Y],npoint[Z],X)*dt;
           A(i,j,k,Y) = (1-dt)*A(i,j,k,Y)+A(npoint[X],npoint[Y],npoint[Z],Y)*dt;
           A(i,j,k,Z) = (1-dt)*A(i,j,k,Z)+A(npoint[X],npoint[Y],npoint[Z],Z)*dt; 
          }
        }         

}

template<int boundary_size>
void construct_exterior_interior_full(const amrex::Geometry geom ,std::vector<std::array<int,3>> &exterior,amrex::Box &interior){
        construct_exterior_interior<X,boundary_size>(geom,exterior,interior);
        construct_exterior_interior<Y,boundary_size>(geom,exterior,interior);
        construct_exterior_interior<Z,boundary_size>(geom,exterior,interior);
        std::sort(exterior.begin(),exterior.end());
        exterior.erase( std::unique( exterior.begin(), exterior.end() ), exterior.end() );
}

typedef void (*F_INTERIOR)(const amrex::Geometry,amrex::Box const&,amrex::Array4<amrex::Real const> const&,amrex::Array4<amrex::Real> const&,double);
typedef void (*F_EXTERIOR)(const amrex::Geometry,amrex::Array4<amrex::Real > const&,std::vector<std::array<int,3>> const&,double);


// Always periodic 
template<F_INTERIOR InteriorF>
void push_ff(const amrex::Geometry geom, amrex::Box const& bx, amrex::Array4<amrex::Real const> const& Source,amrex::Array4<amrex::Real> const& Target,double dt){ 
   auto interior=bx;
   InteriorF(geom,interior,Source,Target,dt);

}
// Periodicity can be set at start
template<F_INTERIOR InteriorF,F_EXTERIOR ExteriorF,int boundary_size=1> 
void push_ff(const amrex::Geometry geom, amrex::Box const& bx, amrex::Array4<amrex::Real const> const& Source,amrex::Array4<amrex::Real> const& Target,double dt){ 
   auto interior=bx;
   std::vector<std::array<int,3>> exterior;

    if(!geom.isAllPeriodic()){
        construct_exterior_interior_full<boundary_size>(geom,exterior,interior);
        ExteriorF(geom,Target,exterior,dt);   
    }

   InteriorF(geom,interior,Source,Target,dt);

}


/* Partial specializations are not possible in c++
 * But here is the code for an arbitrary order iterator, probably no one wants to do more than a couple
 * of orders, so creating them by hand should be ok
 
template<int order,int W_range>
inline void Theta_map(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab &E, amrex::MultiFab &B,double dt ){
   const int l=order/2-1;
   const double alpha=1/(2-pow(2,1/(2*l+1)));
   const double beta=1-2*alpha;
   Theta_map<order-2>(geom,P,E,B,alpha*dt);
   Theta_map<order-2>(geom,P,E,B,beta*dt);
   Theta_map<order-2>(geom,P,E,B,alpha*dt);
}
*/




template<int W_range>
inline void Theta_map1(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab &E, amrex::MultiFab &B,double dt ){

    G_Theta_B(geom,P,E,B,dt);
    G_Theta_E<W_range>(geom,P,E,B,dt);
    G_Theta<Z,W_range>(geom,P,E,B,dt);
    G_Theta<Y,W_range>(geom,P,E,B,dt);
    G_Theta<X,W_range>(geom,P,E,B,dt);

}

template<int W_range>
inline void Theta_map2(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab &E, amrex::MultiFab &B,double dt ){

    G_Theta_E<W_range>(geom,P,E,B,dt/2);
    G_Theta<X,W_range>(geom,P,E,B,dt/2);
    G_Theta<Y,W_range>(geom,P,E,B,dt/2);
    G_Theta<Z,W_range>(geom,P,E,B,dt/2);
    G_Theta_B(geom,P,E,B,dt);
    G_Theta<Z,W_range>(geom,P,E,B,dt/2);
    G_Theta<Y,W_range>(geom,P,E,B,dt/2);
    G_Theta<X,W_range>(geom,P,E,B,dt/2);
    G_Theta_E<W_range>(geom,P,E,B,dt/2);

}

template<int W_range>
inline void Theta_map4(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab &E, amrex::MultiFab &B,double dt ){
   int order=4;
   const int l=order/2-1;
   const double alpha=1/(2-pow(2,1/(2*l+1)));
   const double beta=1-2*alpha;
   Theta_map2<W_range>(geom,P,E,B,alpha*dt);
   Theta_map2<W_range>(geom,P,E,B,beta*dt);
   Theta_map2<W_range>(geom,P,E,B,alpha*dt);
}

#endif
