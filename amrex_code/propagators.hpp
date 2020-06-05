#ifndef PROPAGATOR
#define PROPAGATOR
#include<AMReX_Geometry.H>
#include <array>
#include "w_defs.hpp"
#include "AMReX_Box.H"
#include "AMReX_BoxArray.H"
#include "AMReX_BoxDomain.H"
#include "AMReX_MultiFab.H"
#include "particle_defs.hpp"
#include "amrex_util.hpp"
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


// Global update, they also handle triggering global communication
void G_Theta_E(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab &E, amrex::MultiFab &B,double dt );
void G_Theta_B(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab &E, amrex::MultiFab &B,double dt );

inline int sign(int x) {
    return (x > 0) - (x < 0);
}
inline int mod(int a,int n){
return (a%n+n)%n;
}



// This is for one particle type
// If there are several you need to do this again
// 0 -> x  , 1-> y 2->z
template<int comp,int W_range>
void Theta(CParticles&particles, const amrex::Geometry geom,amrex::Array4<amrex::Real> const& E ,amrex::Array4<amrex::Real> const& B  ,amrex::Box bx ,int ng,double dt){ 
    const auto low = geom.ProbLo();
    const auto Ics = geom.InvCellSize();
    // Remote particle grid lower corner
    for(auto& p : particles){


    const double m= p.rdata(M);
    const double q= p.rdata(Q);
    const  double B_coef = q/m*geom.CellSize(comp);
    const  double E_coef = q*Ics[X]*Ics[Y]*Ics[Z]*geom.CellSize(comp);
        double new_pos=p.pos(comp)+dt*p.rdata(comp+2);
        int num_segments=get_num_segments<comp>(geom,p.pos(comp),new_pos);
       // Never more than two segments!! 
        std::array<std::tuple<amrex::Real,amrex::Real,int>,2> segments;

        get_segment_list<comp>(geom,segments,num_segments, p.pos(comp) , new_pos);
        int coord[3];
        coord[X]=floor((p.pos(X) -low[X])*Ics[X]);
        coord[Y]=floor((p.pos(Y) -low[Y])*Ics[Y]);
        coord[Z]=floor((p.pos(Z) -low[Z])*Ics[Z]);
        

        amrex::Real res_c1=0;
        amrex::Real res_c2=0;


        constexpr int W1_r=W_range*2;
        constexpr int W1_li=-W_range+1;
        constexpr int W1_hi= W_range;
        constexpr int Wp_li= W1_li;
        constexpr int Wp_hi= W1_hi-1;
        constexpr int Wp_r=W1_r-1;

        auto comp_u = (comp+1)%3;
        auto comp_l = (comp+2)%3;
        
            std::array<double,W1_r> comp_uW1={0};
            std::array<double,W1_r> comp_lW1={0};
            std::array<double,W1_r> comp_uWp={0};
            std::array<double,W1_r> comp_lWp={0};




            auto nl=(p.pos(comp_l)-low[comp_l])*Ics[comp_l];
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
            auto nu=(p.pos(comp_u)-low[comp_u])*Ics[comp_u];
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

        for(int seg=0;seg<num_segments;seg++){
            coord[comp] = std::get<2>(segments[seg]);
            auto i_s = std::get<0>(segments[seg]); 
            auto i_e = std::get<1>(segments[seg]);
            

            std::array<double,Wp_r> compI_W12={0};
           
           auto ncs=(i_s-low[comp])*Ics[comp];
           auto nce=(i_e-low[comp])*Ics[comp];
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

                    E(cx,cy,cz,comp)-=E_coef*mul*compI_W12[idc];
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
            if(!geom.isPeriodic(comp)){
                auto res = reflect_boundary<comp>(geom,p.pos(comp)+dt*p.rdata(2+comp));
                p.pos(comp) =res.first;
                p.rdata(comp+2)*=res.second;   
            }
            else{
            p.pos(comp)+=dt*p.rdata(2+comp);
            shift_periodic<comp>(geom,p);
            }
        // B Vel update
        p.rdata( (comp +2 )% 3 +2   )+=B_coef*res_c1;
        p.rdata( (comp+1) % 3 +2  )+=B_coef*res_c2; 

    }
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
    auto coord =get_point_cell(geom,{p.pos(X),p.pos(Y),p.pos(Z)}) ;
       
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



template <int coord>
void G_Theta(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab &E, amrex::MultiFab &B,double dt ){

    E.setBndry(0);
    for (amrex::MFIter mfi(E); mfi.isValid(); ++mfi){
        auto box=E.box(mfi.index());
         
        // Each grid,tile has a their own local particle container
        auto& Part = P.GetParticles(0)[std::make_pair(mfi.index(),mfi.LocalTileIndex())];
        auto&  particles = Part.GetArrayOfStructs();
        amrex::FArrayBox& bfab =B[mfi];
        
        amrex::Array4<amrex::Real> const& B_loc = bfab.array(); 
        amrex::Array4<amrex::Real> const& E_loc = E[mfi].array(); 
        int ng = E.nGrow();
        Theta<coord,WRANGE>(particles,geom,E_loc,B_loc,box,ng,dt);
    }
    
    E.SumBoundary(geom.periodicity());
    P.Redistribute();
    E.FillBoundary(geom.periodicity());
    B.FillBoundary(geom.periodicity());
    

    
}


template<int coord,int boundary_size>
void construct_exterior_interior(const amrex::Geometry geom ,std::vector<std::array<int,3>> &exterior,amrex::Box &interior){
   const auto domain=geom.Domain();
   const auto Lo = domain.loVect();
   const auto Hi = domain.hiVect();
   const auto local_lo=interior.loVect();
   const auto local_hi=interior.hiVect();
   if(!geom.isPeriodic(coord) && boundary_size > 0){
        std::array<int,3> point;
        const int up_c =(coord+1) % 3;
        const int lo_c =(coord+2) % 3;
       if(interior.loVect()[coord]==Lo[coord]){
           interior.growLo(coord,-boundary_size);
            for(int k=local_lo[up_c]; k <= local_hi[up_c];k++){
                for(int j=local_lo[lo_c]; j<=local_hi[lo_c];j++ ){
                    for(int i=Lo[coord]; i < Lo[coord]+boundary_size;i++){
                        point[up_c]=k;
                        point[lo_c]=j;
                        point[coord]=i;
                        exterior.push_back(point);

                    }
                }
            }
       }
        
       if(interior.hiVect()[coord]==Hi[coord]){
           interior.growHi(coord,-boundary_size);
            for(int k=local_lo[up_c]; k <= local_hi[up_c];k++){
                for(int j=local_lo[lo_c]; j<=local_hi[lo_c];j++ ){
                    for(int i=Hi[coord]; i > Hi[coord]-boundary_size;i--){
                        point[up_c]=k;
                        point[lo_c]=j;
                        point[coord]=i;
                        exterior.push_back(point);
                    }
                }
            }
       }
   }     
}

template<int coord>
void MABC(const amrex::Geometry geom,amrex::Array4<amrex::Real> const& A,std::vector<std::array<int,3>> &exterior,double dt){
   const auto domain=geom.Domain();
   const auto Lo = domain.loVect();
   const auto Hi = domain.hiVect();
        for(auto point : exterior){
          if(point[coord] == Lo[coord]){
            const int i = point[X];
            const int j = point[Y];
            const int k = point[Z];
            auto npoint=point;
            npoint[coord]+=1;
           A(i,j,k,X) = (1-dt)*A(i,j,k,X)+A(npoint[X],npoint[Y],npoint[Z],X)*dt;
           A(i,j,k,Y) = (1-dt)*A(i,j,k,Y)+A(npoint[X],npoint[Y],npoint[Z],Y)*dt;
           A(i,j,k,Z) = (1-dt)*A(i,j,k,Z)+A(npoint[X],npoint[Y],npoint[Z],Z)*dt; 
          }
          else if(point[coord] == Hi[coord]){ 
            const int i = point[X];
            const int j = point[Y];
            const int k = point[Z];
            auto npoint=point;
            npoint[coord]-=1;
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

template<int boundary_size,typename InteriorF,typename ExteriorF> 
void push_ff(const amrex::Geometry geom, amrex::Box const& bx, amrex::Array4<amrex::Real const> const& Source,amrex::Array4<amrex::Real> const& Target, InteriorF update_interior,ExteriorF update_exterior,double dt){
   
   auto interior=bx;
   std::vector<std::array<int,3>> exterior;

    if(!geom.isAllPeriodic()){
        construct_exterior_interior_full<boundary_size>(geom,exterior,interior);
        update_exterior(geom,Target,exterior,dt);   
    }

   update_interior(geom,interior,Source,Target,dt);

}

#endif


template<int order>
inline void Theta_map(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab &E, amrex::MultiFab &B,double dt ){
   const int l=order/2-1;
   const double alpha=1/(2-pow(2,1/(2*l+1)));
   const double beta=1-2*alpha;
   Theta_map<order-2>(geom,P,E,B,alpha*dt);
   Theta_map<order-2>(geom,P,E,B,beta*dt);
   Theta_map<order-2>(geom,P,E,B,alpha*dt);
}
template<>
inline void Theta_map<1>(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab &E, amrex::MultiFab &B,double dt ){

    G_Theta_B(geom,P,E,B,dt);
    G_Theta_E(geom,P,E,B,dt);
    G_Theta<Z>(geom,P,E,B,dt);
    G_Theta<Y>(geom,P,E,B,dt);
    G_Theta<X>(geom,P,E,B,dt);

}
template<>
inline void Theta_map<2>(const amrex::Geometry geom,CParticleContainer&P, amrex::MultiFab &E, amrex::MultiFab &B,double dt ){

    G_Theta_E(geom,P,E,B,dt/2);
    G_Theta<X>(geom,P,E,B,dt/2);
    G_Theta<Y>(geom,P,E,B,dt/2);
    G_Theta<Z>(geom,P,E,B,dt/2);
    G_Theta_B(geom,P,E,B,dt);
    G_Theta<Z>(geom,P,E,B,dt/2);
    G_Theta<Y>(geom,P,E,B,dt/2);
    G_Theta<X>(geom,P,E,B,dt/2);
    G_Theta_E(geom,P,E,B,dt/2);

}


