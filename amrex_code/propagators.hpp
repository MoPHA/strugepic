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





void fill_extra_halo(amrex::Geometry geom, amrex::Array4<amrex::Real> const& A,amrex::Box bx,int ng);

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
template<int comp>
void map_from_aux(const amrex::Geometry geom,amrex::Array4<amrex::Real> const& E,amrex::Array4<amrex::Real const> const& E_L ,amrex::Box box_L,amrex::Box box_S,int ng){



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

    // AMReX does not implement halo regions for narrow subgrids
    // I.e. The halo region can not reach over two other subgrids
    // This is ok as you never want very small subgrids 2x2 or 1x1
    // except for when a subgrid shares a periodic boundary with itself
    // 2D slab in 3D space.
    // This is a quick and dirty fix for that
    // So Nx1x1 domains are not working properly at the moment
    
    const auto domainLo = geom.Domain().loVect();
    const auto subgridLo =box_S.loVect();
    const auto domainHi = geom.Domain().hiVect();
    const auto subgridHi = box_S.hiVect();
    const auto auxgridLo = box_L.loVect();
    const auto auxgridHi = box_L.hiVect();

    for(auto coord: {X,Y,Z}){
        auto coord_u = (coord+1)%3;
        auto coord_l = (coord+2)%3;
        auto missing_L =  ng-(domainHi[coord]-domainLo[coord]+1);
        if(domainLo[coord] == subgridLo[coord] && domainHi[coord]==subgridHi[coord] && geom.isPeriodic(coord)
                && missing_L > 0){
          //  std::cout <<"Missing: " <<missing_L << std::endl;
          //  std::cout <<"COORD: " << coord << std::endl;
             for(int u=0; u<=(subgridHi[coord_u]-subgridLo[coord_u]);u++){
            //     std::cout << "U: " << u << std::endl;
                for(int l=0; l<=(subgridHi[coord_l]-subgridLo[coord_l]);l++){
            //     std::cout << "L: " << l << std::endl;
                   for(int c =0; c < missing_L ; c++){
             //          std::cout << "C: "<< c << std::endl;
                        std::array<int,3> ext_coord;
                        ext_coord[coord_u] = u+auxgridLo[coord_u]+(ng);
                        ext_coord[coord_l] = l+auxgridLo[coord_l]+(ng);
                        ext_coord[coord] = c+auxgridHi[coord]-(missing_L-1);
                    std::array<int,3> int_coord;
                    int_coord[coord_l]=ext_coord[coord_l]-shift[coord_l];
                    int_coord[coord_u]=ext_coord[coord_u]-shift[coord_u];
                    int_coord[coord]=ext_coord[coord];
                    
                    int_coord[X]=mod((int_coord[X]-e_low[X]), elx)+e_low[X];
                    int_coord[Y]=mod((int_coord[Y]-e_low[Y]), ely)+e_low[Y];
                    int_coord[Z]=mod((int_coord[Z]-e_low[Z]), elz)+e_low[Z];
            //std::cout << "PATH: " << int_coord[X]<< " "<<int_coord[Y]<< " "<<int_coord[Z]<< " "<<ext_coord[X]<< " "<<ext_coord[Y]<< " "<<ext_coord[Z] <<std::endl;
              
                    E(int_coord[X],int_coord[Y],int_coord[Z],comp)+=E_L(ext_coord[X],ext_coord[Y],ext_coord[Z],comp);
                   } 
                }
             } 
             for(int u=0; u<=(subgridHi[coord_u]-subgridLo[coord_u]);u++){
              //  std::cout << "U: " << u << std::endl;
                for(int l=0; l<=(subgridHi[coord_l]-subgridLo[coord_l]);l++){
              //   std::cout << "L: " << l << std::endl;
                   for(int c =0; c < missing_L ; c++){
                //       std::cout << "C: " << c<< std::endl;
                        std::array<int,3> ext_coord;
                        ext_coord[coord_u] = u+auxgridLo[coord_u]+(ng);
                        ext_coord[coord_l] = l+auxgridLo[coord_l]+(ng);
                        ext_coord[coord] = c+auxgridLo[coord];
                    std::array<int,3> int_coord;
                    int_coord[coord_l]=ext_coord[coord_l]-shift[coord_l];
                    int_coord[coord_u]=ext_coord[coord_u]-shift[coord_u];
                    int_coord[coord]=ext_coord[coord];
                    
                    int_coord[X]=mod((int_coord[X]-e_low[X]), elx)+e_low[X];
                    int_coord[Y]=mod((int_coord[Y]-e_low[Y]), ely)+e_low[Y];
                    int_coord[Z]=mod((int_coord[Z]-e_low[Z]), elz)+e_low[Z];
            //std::cout << "PATH: " << int_coord[X]<< " "<<int_coord[Y]<< " "<<int_coord[Z]<< " "<<ext_coord[X]<< " "<<ext_coord[Y]<< " "<<ext_coord[Z] <<std::endl;
                    E(int_coord[X],int_coord[Y],int_coord[Z],comp)+=E_L(ext_coord[X],ext_coord[Y],ext_coord[Z],comp);
                   } 
                }
             } 
             

        }

    }

}



// This is for one particle type
// If there are several you need to do this again
// 0 -> x  , 1-> y 2->z
template<int comp,int W_range>
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

                    E(cx+shift[X],cy+shift[Y],cz+shift[Z],comp)-=E_coef*mul*compI_W12[idc];
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
        //std::cout << B_loc(0,0,-3,0) <<std::endl;
        //exit(0);
        amrex::Array4<amrex::Real> const& E_L_loc = E_L[mfi].array(); 
        int ng = E_L.nGrow();
         
        Theta<coord,WRANGE>(particles,geom,E_L_loc,B_loc,box_L,box_S,ng,dt);
        fill_extra_halo(geom,B_loc,box_S,E.nGrow());
    }
    E_L.FillBoundary(ggeom.periodicity());
  
    for (amrex::MFIter mfi(E_L); mfi.isValid(); ++mfi){
        auto box_L=mfi.validbox();
        auto box_S=E.box(mfi.index());
        int ng = E_L.nGrow();
        amrex::Array4<amrex::Real > const& E_loc = E.array(mfi);
        amrex::Array4<amrex::Real const> const& E_L_loc = E_L.const_array(mfi); 
        map_from_aux<coord>(geom,E_loc,E_L_loc,box_L,box_S,ng);
        fill_extra_halo(geom,E_loc,box_S,E.nGrow());
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


template<int W_range>
void get_particle_number_density(const amrex::Geometry geom,const amrex::Geometry aux_geom,CParticleContainer&P, amrex::MultiFab &P_dens,amrex::MultiFab &P_dens_aux){
    P_dens.setVal(0);
    P_dens_aux.setVal(0);
    
    int ng=P_dens.nGrow(); 
    for (amrex::MFIter mfi(P_dens_aux); mfi.isValid(); ++mfi){
        auto box_L=mfi.validbox();
        auto box_S=P_dens.box(mfi.index());
    
        // Each grid,tile has a their own local particle container
        auto& Part = P.GetParticles(0)[std::make_pair(mfi.index(),mfi.LocalTileIndex())];
        auto&  particles = Part.GetArrayOfStructs();
        amrex::Array4<amrex::Real> const& P_dens_aux_loc = P_dens_aux[mfi].array(); 

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
        int coord[3];
        coord[X]=floor((p.pos(X) -low[X])*Ics[X]);
        coord[Y]=floor((p.pos(Y) -low[Y])*Ics[Y]);
        coord[Z]=floor((p.pos(Z) -low[Z])*Ics[Z]);
        
            

        auto px=(p.pos(X)-low[X])*Ics[X];
        auto py=(p.pos(Y)-low[Y])*Ics[Y];
        auto pz=(p.pos(Z)-low[Z])*Ics[Z];
        
        constexpr int W1_li=-W_range+1;
        constexpr int W1_hi= W_range;
        constexpr int Wp_li= W1_li;
        constexpr int Wp_hi= W1_hi-1;
        
        for(int i=Wp_li;i<=Wp_hi;i++){
            auto nx=coord[X] +i;
            for(int j=Wp_li;j<=Wp_hi;j++){
                auto ny=coord[Y] +j;
                for(int k=Wp_li;k<=Wp_hi;k++){
                    auto nz=coord[Z] +k;
                       P_dens_aux_loc(nx+shift[X],ny+shift[Y],nz+shift[Z],X)+=
                          Wp(px-nx)*Wp(py-ny)*Wp(pz-nz);
            }
        }
        }

    }

    }
    P_dens_aux.FillBoundary(aux_geom.periodicity());
  
    for (amrex::MFIter mfi(P_dens_aux); mfi.isValid(); ++mfi){
        auto box_L=mfi.validbox();
        auto box_S=P_dens.box(mfi.index());
        amrex::Array4<amrex::Real const > const& P_dens_aux_loc = P_dens_aux.const_array(mfi); 
        amrex::Array4<amrex::Real > const& P_dens_loc = P_dens.array(mfi);
        map_from_aux<0>(geom,P_dens_loc,P_dens_aux_loc,box_L,box_S,ng);
    }

   

}
