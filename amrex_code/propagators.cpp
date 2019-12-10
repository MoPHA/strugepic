#include "AMReX_Box.H"
#include "AMReX_Geometry.H"

void push_B( amrex::Box const& bx,  amrex::Array4<amrex::Real> const& B, amrex::Array4<amrex::Real> const& E,double dt){
   const auto lo = amrex::lbound(bx);
   const auto hi = amrex::ubound(bx);
   for     (int k = lo.z; k <= hi.z; ++k) {
     for   (int j = lo.y; j <= hi.y; ++j) {
       for (int i = lo.x; i <= hi.x; ++i) { 
         B(i,j,k,0) -=dt* (E(i+1,j,k,0)-E(i,j,k,0));
         B(i,j,k,1) -=dt* (E(i,j+1,k,1)-E(i,j,k,1));
         B(i,j,k,2) -=dt* (E(i,j,k+1,2)-E(i,j,k,2));
       }
     }
   }

}

void push_E( amrex::Box const& bx,  amrex::Array4<amrex::Real> const& E, amrex::Array4<amrex::Real> const& B,double dt){
   const auto lo = amrex::lbound(bx);
   const auto hi = amrex::ubound(bx);
   for     (int k = lo.z; k <= hi.z; ++k) {
     for   (int j = lo.y; j <= hi.y; ++j) {
       for (int i = lo.x; i <= hi.x; ++i) { 
         E(i,j,k,0) +=dt* (B(i+1,j,k,0)-B(i,j,k,0));
         E(i,j,k,1) +=dt* (B(i,j+1,k,1)-B(i,j,k,1));
         E(i,j,k,2) +=dt* (B(i,j,k+1,2)-B(i,j,k,2));
       }
     }
   }

}

template<class Particles>
void push_V_E( const Particles&particles, const amrex::Geometry geom,double dt,amrex::Array4<amrex::Real> const& E ){
    for(auto& p : particles){
        auto index =get_point_cell(geom,{p.pos(0),p.pos(1),p.pos(2)}) ;
        // Basis functions are supported over two cells in each direction
        // Need some other particle strucutre to have constant m and q for some
        const int idx_list[4]={-2,1,0,1};
        const double m= p.rdata(0);
        const double q= p.rdata(1);
       
        const double coef = dt*q/m;
        const double inv_dx = 1/geom.CellSize(0);
        const double inv_dy = 1/geom.CellSize(1);
        const double inv_dz = 1/geom.CellSize(2);
        double dvx=0;
        double dvy=0;
        double dvz=0;


        for(auto k: idx_list){
            for(auto j: idx_list){
                for(auto i: idx_list){
                   dvx+= 
                   dvy+= 
                   dvz+= 
                }
            }
        }
        p.rdata(2)+=dvx*coef;
        p.rdata(3)+=dvy*coef;
        p.rdata(4)+=dvz*coef;
    }
}

