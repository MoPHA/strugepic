#include "AMReX_Box.H"

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
void push_V_E( const Particles&particles ){
    for(auto& p : particles){
        p.getreal(3)=1;
    }
}

