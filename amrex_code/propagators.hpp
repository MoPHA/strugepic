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

void push_V_E( CParticles&particles, const amrex::Geometry geom,amrex::Array4<amrex::Real const> const& E ,double dt);

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
template<int comp>
void update_E(const amrex::Geometry geom,amrex::Array4<amrex::Real> const& E,amrex::Array4<amrex::Real const> const& E_L ,amrex::Box box_L,amrex::Box box_S,int ng){

    auto Ll=amrex::lbound(box_L);
    auto Lu=amrex::ubound(box_L);
    const auto el_low=box_L.loVect();
    const auto e_low=box_S.loVect();
    const auto e_high=box_S.hiVect();
    std::array<int,3> shift;
    shift[X]=(el_low[X]+ng)-e_low[X];
    shift[Y]=(el_low[Y]+ng)-e_low[Y];
    shift[Z]=(el_low[Z]+ng)-e_low[Z];


    // Update interior
    for(int k = Ll.z+ng;k<= Lu.z-ng;k++){
    for(int j = Ll.y+ng;j<= Lu.y-ng;j++){
    for(int i = Ll.x+ng;i<= Lu.x-ng;i++){
            E(i-shift[X],j-shift[Y],k-shift[Z],comp)+=E_L(i,j,k,comp);    
            }
        }
    }
    // The exterior corresponds to contributions from particles owned by other processes
    // No T-junctions in the grid are allowed!!
    // Ghost region is not allowed to reach over on whole grid
    
    // We start by getting the lower left corner of the nghost^3 cubes in the corners of the larger array
    // So this is all in the ghost layer of E_L
    
    std::array<int,3> ll={Ll.x-ng,Ll.y-ng,Ll.z-ng};
    // xp,yp,zp,xd,yd,zd
    std::array<std::array<int,6>,8> corners;
    auto xs= ng + (Lu.x-Ll.x+1);
    auto ys= ng + (Lu.y-Ll.y+1);
    auto zs= ng + (Lu.z-Ll.z+1);

    auto elx=e_high[X]-e_low[X] +1;
    auto ely=e_high[Y]-e_low[Y] +1;
    auto elz=e_high[Z]-e_low[Z] +1;

    int n_c=0;
    for(int x=0; x <2;x++){
        for(int y=0; y <2;y++){
            for(int z=0; z <2;z++){
                corners[n_c]={ll[X]+x*xs,ll[Y]+y*ys,ll[Z]+z*zs,sign(x*2-1),sign(y*2-1),sign(z*2-1)};
                n_c++;
            }

        }
    }
    // Iterate over the cubes 
    for(auto p:corners){
    for(int x=0; x <ng;x++){
        for(int y=0; y <ng;y++){
            for(int z=0; z <ng;z++){
               std::array<int,3> ext_coord={p[X]+x,p[Y]+y,p[Z]+z}; 
               std::array<int,3> int_coord={
                   ((ext_coord[X]-p[3+X]*2*ng)-shift[X]),
                   ((ext_coord[Y]-p[3+Y]*2*ng)-shift[Y]),
                   ((ext_coord[Z]-p[3+Z]*2*ng)-shift[Z])
               };
                // Modulo in case we wrap around 
                int_coord[X]=mod((int_coord[X]-e_low[X]), elx)+e_low[X];
                int_coord[Y]=mod((int_coord[Y]-e_low[Y]), ely)+e_low[Y];
                int_coord[Z]=mod((int_coord[Z]-e_low[Z]), elz)+e_low[Z];
                E(int_coord[X],int_coord[Y],int_coord[Z],comp)+=E_L(ext_coord[X],ext_coord[Y],ext_coord[Z],comp);
            }
        }
    }
    }


    // Construct a list of the lower left corner of the  6 faces 
    // There probably is some better way
    std::array<std::array<int,6>,6> faces;
    faces[0]={Ll.x-ng,Ll.y+ng,Ll.z+ng,-1,0,0};
    faces[1]={faces[0][X]+xs,faces[0][Y],faces[0][Z],1,0,0};
    faces[2]={Ll.x+ng,Ll.y-ng,Ll.z+ng,0,-1,0};
    faces[3]={faces[2][X],faces[2][Y]+ys,faces[0][Z],0,1,0};
    faces[4]={Ll.x+ng,Ll.y+ng,Ll.z-ng,0,0,-1};
    faces[5]={faces[4][X],faces[4][Y],faces[4][Z]+zs,0,0,1};
    
    // Loop limits, again this could be generated
    std::array<std::array<int,3>,6> loop_lims;
    loop_lims[0]={ng,ely,elz};
    loop_lims[1]={ng,ely,elz};
    loop_lims[2]={elx,ng,elz};
    loop_lims[3]={elx,ng,elz};
    loop_lims[4]={elx,ely,ng};
    loop_lims[5]={elx,ely,ng};


    for( int f=0; f<6; f++  ){
    for(int x=0; x <loop_lims[f][X];x++){
        for(int y=0; y <loop_lims[f][Y];y++){
            for(int z=0; z <loop_lims[f][Z];z++){
               std::array<int,3> ext_coord={faces[f][X]+x,faces[f][Y]+y,faces[f][Z]+z}; 
               std::array<int,3> int_coord={
                            (ext_coord[X]-2*ng*faces[f][X+3]-shift[X]),                           
                            (ext_coord[Y]-2*ng*faces[f][Y+3]-shift[Y]),                            
                            (ext_coord[Z]-2*ng*faces[f][Z+3]-shift[Z])                           
                        };

                // Modulo in case we wrap around 
                            int_coord[X]=mod((int_coord[X]-e_low[X]), elx)+e_low[X];
                            int_coord[Y]=mod((int_coord[Y]-e_low[Y]), ely)+e_low[Y];
                            int_coord[Z]=mod((int_coord[Z]-e_low[Z]), elz)+e_low[Z];
                E(int_coord[X],int_coord[Y],int_coord[Z],comp)+=E_L(ext_coord[X],ext_coord[Y],ext_coord[Z],comp);
                }
            }
        }
    }

//Iterate over the edges of the cube
    // Starting points
    std::array<std::array<int,6>,12> edges;
    edges[0]={Ll.x+ng,Ll.y-ng,Ll.z-ng,0,-1,-1};
    edges[1]={edges[0][X],edges[0][Y],edges[0][Z]+zs,0,-1,1};
    edges[2]={edges[0][X],edges[0][Y]+ys,edges[0][Z],0,1,-1};
    edges[3]={edges[0][X],edges[0][Y]+ys,edges[0][Z]+zs,0,1,1};
    
    edges[4]={Ll.x-ng,Ll.y+ng,Ll.z-ng,-1,0,-1};
    edges[5]={edges[4][X],edges[4][Y],edges[4][Z]+zs,-1,0,1};
    edges[6]={edges[4][X]+xs,edges[4][Y],edges[4][Z],1,0,-1};
    edges[7]={edges[4][X]+xs,edges[4][Y],edges[4][Z]+zs,1,0,1};
    
    edges[8]={Ll.x-ng,Ll.y-ng,Ll.z+ng,-1,-1,0};
    edges[9]={edges[8][X]+xs,edges[8][Y],edges[8][Z],1,-1,0};
    edges[10]={edges[8][X],edges[8][Y]+ys,edges[8][Z],-1,1,0};
    edges[11]={edges[8][X]+xs,edges[8][Y]+ys,edges[8][Z],1,1,0};
    // Loop limits
    std::array<std::array<int,3>,12> ed_lims;

    ed_lims[0]={elx,ng,ng};
    ed_lims[1]={elx,ng,ng};
    ed_lims[2]={elx,ng,ng};
    ed_lims[3]={elx,ng,ng};

    ed_lims[4]={ng,ely,ng};
    ed_lims[5]={ng,ely,ng};
    ed_lims[6]={ng,ely,ng};
    ed_lims[7]={ng,ely,ng};

    ed_lims[8]={ng,ng,elz};
    ed_lims[9]={ng,ng,elz};
    ed_lims[10]={ng,ng,elz};
    ed_lims[11]={ng,ng,elz};


    for( int e=0; e<12; e++  ){
    for(int x=0; x <ed_lims[e][X];x++){
        for(int y=0; y <ed_lims[e][Y];y++){
            for(int z=0; z <ed_lims[e][Z];z++){
               std::array<int,3> ext_coord={edges[e][X]+x,edges[e][Y]+y,edges[e][Z]+z}; 
               std::array<int,3> int_coord={
                            (ext_coord[X]-2*ng*edges[e][X+3]-shift[X]),                           
                            (ext_coord[Y]-2*ng*edges[e][Y+3]-shift[Y]),                            
                            (ext_coord[Z]-2*ng*edges[e][Z+3]-shift[Z])                           
                        };

                // Modulo in case we wrap around 
                            int_coord[X]=mod((int_coord[X]-e_low[X]), elx)+e_low[X];
                            int_coord[Y]=mod((int_coord[Y]-e_low[Y]), ely)+e_low[Y];
                            int_coord[Z]=mod((int_coord[Z]-e_low[Z]), elz)+e_low[Z];
                E(int_coord[X],int_coord[Y],int_coord[Z],comp)+=E_L(ext_coord[X],ext_coord[Y],ext_coord[Z],comp);
                }
            }
        }
    }


}



// This is for one particle type
// If there are several you need to do this again
// 0 -> x  , 1-> y 2->z
// Local and neighbour particle list are not the same data structure
template<int comp>
void Theta(CParticles&particles, const amrex::Geometry geom,amrex::Array4<amrex::Real> const& E ,amrex::Array4<amrex::Real> const& B  ,amrex::Box box_L,amrex::Box box_S,int ng,double dt){ 
    const auto low = geom.ProbLo();
    const auto Ics = geom.InvCellSize();
    // Remote particle grid lower corner
     const auto el_low=box_L.loVect();
     const auto e_low=box_S.loVect();
      std::array<int,3> shift;
      shift[X]=(el_low[X]+ng)-e_low[X];
      shift[Y]=(el_low[Y]+ng)-e_low[Y];
      shift[Z]=(el_low[Z]+ng)-e_low[Z];

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

        auto comp_u = (comp+1)%3;
        auto comp_l = (comp+2)%3;
        
            std::array<double,4> comp_uW1={0,0,0,0};
            std::array<double,4> comp_lW1={0,0,0,0};
            std::array<double,4> comp_uW12={0,0,0,0};
            std::array<double,4> comp_lW12={0,0,0,0};




            auto nl=(p.pos(comp_l)-low[comp_l])*Ics[comp_l];
            for(int l=0; l<4;l++){
                auto cl=coord[comp_l]+(l-1); 
                comp_lW1[l]=W1(nl-cl);

            }
            for(int l=0; l<3;l++){
                auto cl=coord[comp_l]+(l-1); 
                comp_lW12[l]=Wp(nl-cl);

            }
            auto nu=(p.pos(comp_u)-low[comp_u])*Ics[comp_u];
            for(int  u=0; u<4;u++){
                auto cu=coord[comp_u]+(u-1); 
                comp_uW1[u]=W1(nu-cu);
            }
            for(int  u=0; u<3;u++){
                auto cu=coord[comp_u]+(u-1); 
                comp_uW12[u]=Wp(nu-cu);
            }

        for(int seg=0;seg<num_segments;seg++){
            coord[comp] = std::get<2>(segments[seg]);
            auto i_s = std::get<0>(segments[seg]); 
            auto i_e = std::get<1>(segments[seg]);
            

            std::array<double,4> compI_W12={0,0,0};
           
           auto ncs=(i_s-low[comp])*Ics[comp];
           auto nce=(i_e-low[comp])*Ics[comp];
            for(int c=0;c<3;c++){
                auto cc=coord[comp]+(c-1);
                compI_W12[c]=I_Wp(ncs- cc ,nce- cc);

            }



            for(int l=0;l<4; l++ ){
                for(int u=0;u<4;u++){
                double mul=comp_lW1[l]*comp_uW1[u];
                double mulu12l1=comp_uW12[u]*comp_lW1[l];
                double mull12u1=comp_lW12[l]*comp_uW1[u];
                for(int c=0;c<3;c++){
                    int cx,cy,cz;

                    if(comp==0){
                     cx=coord[X]+(c-1);
                     cy=coord[Y]+(u-1);
                     cz=coord[Z]+(l-1);
                    }else if(comp==1){
                     cx=coord[X]+(l-1);
                     cy=coord[Y]+(c-1);
                     cz=coord[Z]+(u-1);

                    }
                    else{
                     cx=coord[X]+(u-1);
                     cy=coord[Y]+(l-1);
                     cz=coord[Z]+(c-1);
                    }

                    E(cx+shift[X],cy+shift[Y],cz+shift[Z],comp)-=E_coef*mul*compI_W12[c];
                    res_c1+=B(cx,cy,cz,comp_u)*mull12u1*compI_W12[c];
                    res_c2-=B(cx,cy,cz,comp_l)*compI_W12[c]*mulu12l1;
                }
            }
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





template <int coord>
void G_Theta(const amrex::Geometry geom,const amrex::Geometry ggeom,CParticleContainer&P, amrex::MultiFab &E,amrex::MultiFab &E_L, amrex::MultiFab &B,double dt ){

    E_L.setVal(0) ;   
    for (amrex::MFIter mfi(E_L); mfi.isValid(); ++mfi){
        auto box_L=mfi.validbox();
        auto box_S=E.box(mfi.index());
    
        // Each grid,tile has a their own local particle container
        auto& Part = P.GetParticles(0)[std::make_pair(mfi.index(),mfi.LocalTileIndex())];
        auto&  particles = Part.GetArrayOfStructs();
        amrex::FArrayBox& bfab =B[mfi];
        
        amrex::Array4<amrex::Real> const& B_loc = bfab.array(); 
        amrex::Array4<amrex::Real> const& E_L_loc = E_L[mfi].array(); 
        int ng = E_L.nGrow();
        Theta<coord>(particles,geom,E_L_loc,B_loc,box_L,box_S,ng,dt);
    }
    E_L.FillBoundary(ggeom.periodicity());
  
    for (amrex::MFIter mfi(E_L); mfi.isValid(); ++mfi){
        auto box_L=mfi.validbox();
        auto box_S=E.box(mfi.index());
        int ng = E_L.nGrow();
        amrex::Array4<amrex::Real > const& E_loc = E.array(mfi);
        amrex::Array4<amrex::Real const> const& E_L_loc = E_L.const_array(mfi); 

        update_E<coord>(geom,E_loc,E_L_loc,box_L,box_S,ng);
    }



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
inline void Theta_map(const amrex::Geometry geom,const amrex::Geometry ggeom,CParticleContainer&P, amrex::MultiFab &E,amrex::MultiFab &E_L, amrex::MultiFab &B,double dt ){
   const int l=order/2-1;
   const double alpha=1/(2-pow(2,1/(2*l+1)));
   const double beta=1-2*alpha;
   Theta_map<order-2>(geom,ggeom,P,E,E_L,B,alpha*dt);
   Theta_map<order-2>(geom,ggeom,P,E,E_L,B,beta*dt);
   Theta_map<order-2>(geom,ggeom,P,E,E_L,B,alpha*dt);
}
template<>
inline void Theta_map<1>(const amrex::Geometry geom,const amrex::Geometry ggeom,CParticleContainer&P, amrex::MultiFab &E,amrex::MultiFab &E_L, amrex::MultiFab &B,double dt ){

    G_Theta_B(geom,P,E,B,dt);
    G_Theta_E(geom,P,E,B,dt);
    G_Theta<Z>(geom,ggeom,P,E,E_L,B,dt);
    G_Theta<Y>(geom,ggeom,P,E,E_L,B,dt);
    G_Theta<X>(geom,ggeom,P,E,E_L,B,dt);

}
template<>
inline void Theta_map<2>(const amrex::Geometry geom,const amrex::Geometry ggeom,CParticleContainer&P, amrex::MultiFab &E,amrex::MultiFab &E_L, amrex::MultiFab &B,double dt ){

    G_Theta_E(geom,P,E,B,dt/2);
    G_Theta<X>(geom,ggeom,P,E,E_L,B,dt/2);
    G_Theta<Y>(geom,ggeom,P,E,E_L,B,dt/2);
    G_Theta<Z>(geom,ggeom,P,E,E_L,B,dt/2);
    G_Theta_B(geom,P,E,B,dt);
    G_Theta<Z>(geom,ggeom,P,E,E_L,B,dt/2);
    G_Theta<Y>(geom,ggeom,P,E,E_L,B,dt/2);
    G_Theta<X>(geom,ggeom,P,E,E_L,B,dt/2);
    G_Theta_E(geom,P,E,B,dt/2);

}


