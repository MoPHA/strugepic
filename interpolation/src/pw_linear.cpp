#include<math.h>
#define REAL double



REAL W1(REAL x){
    const auto fx=fabs(x);
    if(fx >=1 ){
            return 0;
    }
    else{
        return 1-fx;
    }
}

REAL W(REAL x,REAL y, REAL z){
    return W1(x)*W1(y)*W1(z);
}

REAL W1d(REAL x){
    if(x <1 && x >0){
        return -1;
    }
    else if( x > -1 &&  x< 0){
        return 1;
    }
    else{
        return 0;
    }

}

REAL IH_W1(REAL x){
    if(x > 1){
        return 1;
    }
    else if(x < -1){
        return 0;
    }
    else{
        return -(x*fabs(x)-2*x-1)*0.5;
    }
}

REAL I_W1(REAL a,REAL b ){
    return IH_W1(b)-IH_W1(a); 
}

REAL Wp(REAL x){
    if(x <1 && x >0){
        return 1;
    }
    else{
        return 0;
    }
}

REAL IH_wp(REAL x){
    if(x < 0){
    return 0;
    }
    else if(x > 1){
     return 1;       
    }
    else{
    return x;
    }
}

REAL I_Wp(REAL a,REAL b){
    return IH_wp(b)-IH_wp(a);
}

