#ifndef PROPAGATOR
#define PROPAGATOR
#include<AMReX_Geometry.H>
#include <array>
#include "../include/interpolation.hpp"
#include "AMReX_Box.H"
#include "AMReX_BoxArray.H"
#include "AMReX_BoxDomain.H"
#include "AMReX_MultiFab.H"
#include "particle_defs.hpp"
#include "amrex_util.hpp"



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
    const int idx_list[4]={-1,0,1,2};
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
    const  double E_coef = q*geom.InvCellSize(X)*geom.InvCellSize(Y)*geom.InvCellSize(Z)*geom.CellSize(comp);
        const auto p_segments = get_segment_list<comp>(geom, p.pos(comp) , p.pos(comp)+dt*p.rdata(comp+2));
        auto coord =get_point_cell(geom,{p.pos(X),p.pos(Y),p.pos(Z)}) ;
        amrex::Real res_c1=0;
        amrex::Real res_c2=0;

        auto comp_u = (comp+1)%3;
        auto comp_l = (comp+2)%3;
        
            std::array<double,4> comp_uW1={0,0,0,0};
            std::array<double,4> comp_lW1={0,0,0,0};
            std::array<double,4> comp_uW12={0,0,0,0};
            std::array<double,4> comp_lW12={0,0,0,0};



            using namespace strugepic;

            auto nl=(p.pos(comp_l)-low[comp_l])*Ics[comp_l];
            for(auto l: idx_list){
                auto cl=coord[comp_l]+l; 
                comp_lW1[l+1]=W1(nl-cl);
                comp_lW12[l+1]=W12(nl-cl);

            }
            auto nu=(p.pos(comp_u)-low[comp_u])*Ics[comp_u];
            for(auto u: idx_list){
                auto cu=coord[comp_u]+u; 
                comp_uW1[u+1]=W1(nu-cu);
                comp_uW12[u+1]=W12(nu-cu);
            }

        for(auto seg: p_segments){
            coord[comp] = std::get<2>(seg);
            auto i_s = std::get<0>(seg); 
            auto i_e = std::get<1>(seg);
            int cl[3];
            

            std::array<double,4> compI_W12={0,0,0,0};
           
           auto ncs=(i_s-low[comp])*Ics[comp];
           auto nce=(i_e-low[comp])*Ics[comp];
            for(auto c: idx_list){
                auto cc=coord[comp]+c;
                compI_W12[c+1]=I_W12(ncs- cc ,nce- cc);

            }

            std::array<int,3> idx;
            for(auto k: idx_list){
                        cl[Z] = coord[Z]+k;
                        idx[Z]=k;
                        
                for(auto j: idx_list){
                        cl[Y] = coord[Y]+j;
                        idx[Y]=j;
                    for(auto i: idx_list){
                        cl[X] = coord[X]+i;
                        idx[X]=i;
                            E(cl[X]+shift[X],cl[Y]+shift[Y],cl[Z]+shift[Z],comp)+=
                                -E_coef
                                *compI_W12[idx[comp]+1]
                                *comp_uW1[idx[comp_u]+1]
                                *comp_lW1[idx[comp_l]+1];

                    }
                }
            }

            for(auto k: idx_list){
                        cl[Z] = coord[Z]+k;
                        idx[Z]=k;
                        
                for(auto j: idx_list){
                        cl[Y] = coord[Y]+j;
                        idx[Y]=j;
                    for(auto i: idx_list){
                        cl[X] = coord[X]+i;
                        idx[X]=i;
                        res_c1+=B(cl[X],cl[Y],cl[Z],comp_u)
                        *compI_W12[idx[comp]+1]
                        *comp_uW1[idx[comp_u]+1]
                        *comp_lW12[idx[comp_l]+1];
                        
                        res_c2-=B(cl[X],cl[Y],cl[Z],comp_l)
                        *compI_W12[idx[comp]+1]
                        *comp_uW12[idx[comp_u]+1]
                        *comp_lW1[idx[comp_l]+1];
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



template<int comp>
void push_pos_pos(CParticleContainer&p_container ,CParticles&local_particles,const amrex::Geometry geom,double dt){
        for(auto &p: local_particles){
            if(!geom.isPeriodic(comp)){
                auto res = reflect_boundary<comp>(geom,p.pos(comp)+dt*p.rdata(2+comp));
                p.pos(comp) =res.first;
                p.rdata(comp+2)*=res.second;   
            }
            else{
            p.pos(comp)+=dt*p.rdata(2+comp);
            shift_periodic<comp>(geom,p);
            }
        }
}


template<int comp>
void push_B_pos(CParticles&particles, const amrex::Geometry geom, const amrex::Array4<amrex::Real> & B ,double dt){ 
    const int idx_list[4]={-1,0,1,2};
    const auto low = geom.ProbLo();
    const auto Ics = geom.InvCellSize();
    // A scaling factor is needed here if dx,dy,dz are not equal!

    for(auto& p : particles){
    const double m= p.rdata(M);
    const double q= p.rdata(Q);
    const  double coef = q/m*geom.CellSize(comp);
        const auto p_segments = get_segment_list<comp>(geom, p.pos(comp) , p.pos(comp)+dt*p.rdata(comp+2));
        auto coord =get_point_cell(geom,{p.pos(X),p.pos(Y),p.pos(Z)}) ;
        amrex::Real res_c1=0;
        amrex::Real res_c2=0;

        int cl[3];

        auto comp_u = (comp+1)%3;
        auto comp_l = (comp+2)%3;
        
            std::array<double,4> comp_uW1={0,0,0,0};
            std::array<double,4> comp_lW1={0,0,0,0};
            std::array<double,4> comp_uW12={0,0,0,0};
            std::array<double,4> comp_lW12={0,0,0,0};



            using namespace strugepic;

            auto nl=(p.pos(comp_l)-low[comp_l])*Ics[comp_l];
            for(auto l: idx_list){
                auto cl=coord[comp_l]+l; 
                comp_lW1[l+1]=W1(nl-cl);
                comp_lW12[l+1]=W12(nl-cl);

            }
            auto nu=(p.pos(comp_u)-low[comp_u])*Ics[comp_u];
            for(auto u: idx_list){
                auto cu=coord[comp_u]+u; 
                comp_uW1[u+1]=W1(nu-cu);
                comp_uW12[u+1]=W12(nu-cu);
            }

        for(auto seg: p_segments){
            coord[comp] = std::get<2>(seg);
            auto i_s = std::get<0>(seg); 
            auto i_e = std::get<1>(seg);
            
            std::array<double,4> compI_W12={0,0,0,0};
           auto ncs=(i_s-low[comp])*Ics[comp];
           auto nce=(i_e-low[comp])*Ics[comp];
            for(auto c: idx_list){
                auto cc=coord[comp]+c;
                compI_W12[c+1]=I_W12(ncs- cc ,nce- cc);

            }
            
                std::array<int,3> idx;
            for(auto k: idx_list){
                        cl[Z] = coord[Z]+k;
                        idx[Z]=k;
                        
                for(auto j: idx_list){
                        cl[Y] = coord[Y]+j;
                        idx[Y]=j;
                    for(auto i: idx_list){
                        cl[X] = coord[X]+i;
                        idx[X]=i;
        
                        
                        res_c1+=B(cl[X],cl[Y],cl[Z],comp_u)
                        *compI_W12[idx[comp]+1]
                        *comp_uW1[idx[comp_u]+1]
                        *comp_lW12[idx[comp_l]+1];
                        
                        res_c2-=B(cl[X],cl[Y],cl[Z],comp_l)
                        *compI_W12[idx[comp]+1]
                        *comp_uW12[idx[comp_u]+1]
                        *comp_lW1[idx[comp_l]+1];

                    } // over i
                } // over j
            } // over k  
        } // over segments

        p.rdata( (comp +2 )% 3 +2   )+=coef*res_c1;
        p.rdata( (comp+1) % 3 +2  )+=coef*res_c2; 

    } // over particles
} // function end




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

#endif

