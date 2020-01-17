
// 8th order piecewise polynominal coefficients, from highest to lowest, comments indicate the range
// Function is zero outside the range
static const double P8_coeffs[36]={
        15.0/1024,15.0/128,49.0/128,21.0/32,35.0/64,0,0,1,1,  // -2 -> -1
        -15.0/1024,15.0/128,7.0/16,21.0/32,175.0/256,0,-105.0/128,0,337.0/512, // -1 -> 0
        -15.0/1024,-15.0/128,7.0/16,-21.0/32,175.0/256,0,-105.0/128,0,337.0/512, // 0 -> 1
        15.0/1024,-15.0/128,49.0/128,-21.0/32,35.0/64,0,0,-1,1 //1->2
        };
// Derivative for previous
static const double P8d_coeffs[32]={
         15.0/128, 105.0/128, 147.0/64,105.0/32,35.0/16,0,0,1, // -2 -> -1 
         -15.0/128,105.0/128,21.0/8,105.0/32,175.0/64,0,-105.0/64,0,  // -1 -> 0
         -15.0/128,-105.0/128,21.0/8,-105.0/32,175.0/64,0,-105.0/64,0,  // 0  -> 1
         15.0/128, -105.0/128, 147.0/64,-105.0/32 ,35.0/16,0,0,-1 // 1 -> 2 
         };

// Integral coefficients (with extra constant so continuous)
static const double P8I_coeffs[40]={
        5.0/3072,15.0/1024,7.0/128,7.0/64,7.0/64,0,0,1.0/2,1,7.0/12, // -2 -> 1
        -5./3072,15./1024,1./16,7./64,35./256,0,-35./128,0,337./512,1.0/2 , // -1->0
        -5./3072,-15./1024,1./16,-7./64,35./256,0,-35./128,0,337./512,1.0/2 , // 0->1
        5.0/3072,-15.0/1024,7.0/128,-7.0/64,7.0/64,0,0,-1.0/2,1,5.0/12 // 1 -> 2
        };

// Polynominal coefficients for different orders low -> High
static const double* P_COEFS[9]={0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,P8_coeffs};
// Coefficients for the derivatives low -> high
static const double* Pd_COEFS[8]={0x0,0x0,0x0,0x0,0x0,0x0,0x0,P8d_coeffs};
// Coefficients for the integrals low -> high
static const double* PI_COEFS[10]={0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,P8I_coeffs};


namespace strugepicInternal 
{
    
    // Function to evaluate  piecewise compact polynominal over -2,2 consisting of 4 parts  
    template <class T,int Pdeg>
    inline T P4_eval(T x,const T* coeffs){
            int start_idx = ( (int )(x)+2)*(Pdeg+1);
            T sum =coeffs[start_idx];
            for(int i=1; i<(Pdeg+1);i++){
                sum = coeffs[start_idx +i]+x*sum;
            }
            return sum;
        }
    

    template <class T,int Pdeg>
        T W1(T x){

        if( x >= 2 ||  x <= -2){
            return 0;
        }
        else{ 
            return P4_eval<T,Pdeg>(x,(T *)P_COEFS[Pdeg]);
        }
        }

    template <class T, int Pdeg>
        T W1d(T x){

        if( x >= 2 ||  x <= -2){
            return 0;
        }
        else{ 
            return P4_eval<T,Pdeg-1>(x,(T *)Pd_COEFS[Pdeg-1]);
            }
        }
    
    template <class T,int Pdeg>
        T W12(T x){
            if(x >= 2 || x <= -1   ){
                return 0;
            }else{
                return -1*(W1d<T,Pdeg>(x)+W1d<T,Pdeg>(x+1)+W1d<T,Pdeg>(x+2));
            }
        }

    
    template <class T, int Pdeg>
        inline T IH_W1(T q){
            if(q >= 2){
                return 1;
            }else if(q < -2){

                return 0;
            }
            else{
                return P4_eval<T,Pdeg+1>(q,(T *)PI_COEFS[Pdeg+1]);
            }
        }

    template <class T, int Pdeg>
        T I_W1(T a, T b){
            return IH_W1<T,Pdeg>(b)-IH_W1<T,Pdeg>(a);
        }


    // Integral function of W12
    template <class T,int Pdeg>
      inline  T IH_W12(T q){
            if(q > 2 ){
                return 1;
            }
            else if(q< -1){
                return 0;
            }
            else{
                return -1*(W1<T,Pdeg>(q)+W1<T,Pdeg>(q+1)+W1<T,Pdeg>(q+2))+1;
            }
        }



    // Integrate over a range.
    template <class T,int Pdeg>
        T I_W12(T a, T b ){
          return IH_W12<T,Pdeg>(b) - IH_W12<T,Pdeg>(a); 
        }

}
