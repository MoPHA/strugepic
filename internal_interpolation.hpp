

namespace strugepicInternal 
{
    
    template <class T>
    T W1_8 ( T x ){
        static const T coeffs[36]={
        15.0/1024,15.0/128,49.0/128,21.0/32,35.0/64,0,0,1,1,  // -2 -> -1
        -15.0/1024,15.0/128,7.0/16,21.0/32,175.0/256,0,-105.0/128,0,337.0/512, // -1 -> 0
        -15.0/1024,-15.0/128,7.0/16,-21.0/32,175.0/256,0,-105.0/128,0,337.0/512, // 0 -> 1
        15.0/1024,-15.0/128,49.0/128,-21.0/32,35.0/64,0,0,-1,1 //1->2
        };

        
        if( x >= 2 ||  x <= -2){
            return 0;
        }
        else{ 
            int start_idx = ( (int )(x+2))*9;
            T sum =coeffs[start_idx];
            for(int i=1; i<9;i++){
                sum = coeffs[start_idx +i]+x*sum;
            }
            return sum;
        }
    
    }



    

    template <class T>
    T W1d_8 (T q){
        static const T coeffs[32]={
         15.0/128, 105.0/128, 147.0/64,105.0/32,35.0/16,0,0,1, // -2 -> -1 
         -15.0/128,105.0/128,21.0/8,105.0/32,175.0/64,0,-105.0/64,0,  // -1 -> 0
         -15.0/128,-105.0/128,21.0/8,-105.0/32,175.0/64,0,-105.0/64,0,  // 0  -> 1
         15.0/128, -105.0/128, 147.0/64,-105.0/32 ,35.0/16,0,0,-1 // 2 -> 1 
        };

        if( q >= 2 ||  q <= -2){
            return 0;
        }
        else{
            int start_idx = ( (int )(q+2))*8;
            T sum =coeffs[start_idx];
            for(int i=1; i<8;i++){
                sum = coeffs[start_idx +i]+q*sum;
            }
            return sum;
        }

    }

    template <class T>
    T W12_8 (T q){
        if(q >= 2 || q <= -1 ){
            return 0;
        }
        else{
            return W1d_8<T>(q)+W1d_8<T>(q+1)+W1d_8<T>(q+2);
        }

        
    }


 
 
}
