#ifndef POLYUTIL
#define POLYUTIL

#ifdef GPUBUILD
    #define DEVICE_QUALIFIER __host__ __device__  
#else
    #define DEVICE_QUALIFIER
#endif
// Some quick utility functions 
// To evaluate piecewise polynomials in a compact region
// with integer limits and integer bins
// Polynomial is also normalized over the range

// Only meant for compiling the default interpolation function
// and not exposed anywhere else, that's why it's not placed in the 
// top include directory


 template <class T,int Poly_degree,int lower_limit>
   DEVICE_QUALIFIER   inline T P_eval(T x,const T* coeffs){
            int start_idx = ( floor(x)-lower_limit)*(Poly_degree+1);
            T sum =coeffs[start_idx];
            for(int i=1; i<(Poly_degree+1);i++){
                sum = coeffs[start_idx +i]+x*sum;
            }
            return sum;
        }

 template <class T,int Poly_degree,int lower_limit,int upper_limit>
   DEVICE_QUALIFIER  inline T P(T x,const T* coeffs){

        if( x >= upper_limit ||  x <= lower_limit){
            return 0;
        }
        else{ 
            return P_eval<T,Poly_degree,lower_limit>(x,coeffs);
        }
        }

  template <class T, int Poly_degree,int lower_limit, int upper_limit>
   DEVICE_QUALIFIER   inline T IP(T q,const T* coeffs){
            if(q >= upper_limit){
                return 1;
            }else if(q < lower_limit){

                return 0;
            }
            else{
                return P_eval<T,Poly_degree+1,lower_limit>(q,coeffs);
            }
        }

#endif
