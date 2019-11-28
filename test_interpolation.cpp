#include<iostream>
#include"interpolation.hpp"
#include<iomanip>
#include <vector>
#include <fstream>
#include<cmath>



// Returns 0 if All is ok
int compare_res(double ref_val, double calc_val){
        double max_diff=1e-9;
        int res; 
        res= !(fabs(calc_val-ref_val) < max_diff);
        if(res !=0){
            std::cerr << "Error is to large, Value is: " <<calc_val << " Should be:"<< ref_val << std::endl;  
        }
        return res;
};


template< int NumArgs>
double ev_f(void * test_func, std::vector<double>* coords,int i){

double (*tf1)(double);
double (*tf2)(double,double);
double (*tf3)(double,double,double);

if(NumArgs ==1){
    tf1 = (*(double(**)(double)) (&test_func)); 
    return tf1(coords[0][i]);
}
else if(NumArgs ==2){
    tf2 = (*(double(**)(double,double)) (&test_func));
    return tf2(coords[0][i],coords[1][i]);
}
else if(NumArgs ==3){
    tf3 = (*(double(**)(double,double,double)) (&test_func));
    return tf3(coords[0][i],coords[1][i],coords[2][i]);
}
else{

    return -2000000;
}

}


template<int NumArgs >
int run_test(std::string test_name,void * test_func,std::string *arg_data_files,std::string ref_data_file  ){

    std::cout << "TEST: " << test_name;
    std::cerr << "TEST: " << test_name << std::endl;
    
    double num=0.0;
    std::vector<double> coord_vals[NumArgs];
    std::vector<double> ref_vals;
    std::ifstream coordfiles[NumArgs];
    
    for(int i=0; i< NumArgs ; i++){
        coordfiles[i].open(arg_data_files[i], std::ios::in);
        if(!coordfiles[i].is_open()){
            return 2;
        }
    }


    std::ifstream reffile(ref_data_file, std::ios::in);
    if(!reffile.is_open()){
        return 2;
    }



    for(int i=0; i<NumArgs; i++){
        while (coordfiles[i] >> num) {
            coord_vals[i].push_back(num);
        }
    } 

    while (reffile >> num) {
            ref_vals.push_back(num);
    }

    int cr;

    
    for(int i=0; i< (int) ref_vals.size()-1;i++){
        
        double f_res=ev_f<NumArgs>(test_func,coord_vals,i);

        cr=compare_res(ref_vals[i],f_res);
        if(cr !=0){
            std::cout << " FAIL"  << std::endl;
            return 1;
        }
    }

    
    std::cout << " PASS" << std::endl;
    return 0 ;
}



int main(void){
    std::cout << "Running tests:" << std::endl;

    std::string coordf[3]={"data_x.out","data_y.out","data_z.out"};

    auto f1= strugepic::W;
    void *tf1 = *(void**)(&f1);
    run_test<3>("W(x,y,z)",tf1,coordf,"data_W.out");

    auto f2= strugepic::W1d;
    void *tf2 = *(void**)(&f2);
    run_test<1>("W1d(x)",tf2,coordf,"./data_W1d.out");

    auto f3= strugepic::W12;
    void *tf3 = *(void**)(&f3);
    run_test<1>("W12(x)",tf3,coordf,"./data_W12.out");
    
    auto f4= strugepic::I_W12;
    void *tf4 = *(void**)(&f4);
    run_test<2>("IW12",tf4,coordf,"./data_IW12.out");

}
