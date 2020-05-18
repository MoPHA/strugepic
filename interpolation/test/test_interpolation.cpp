#include<iostream>
#include"w_defs.hpp"
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
            std::cerr << std::setprecision(16) << "Error is to large, Value is: " <<calc_val << " Should be:"<< ref_val << std::endl;  
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
            std::cerr << "Could not open file: " << arg_data_files[i] << std::endl;
            std::cout << " FAIL"  << std::endl;
            return 2;
        }
    }


    std::ifstream reffile(ref_data_file, std::ios::in);
    if(!reffile.is_open()){
        std::cerr << "Could not open file: " << ref_data_file << std::endl;
        std::cout << " FAIL"  << std::endl;
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

    auto f1= W;
    void *tf1 = *(void**)(&f1);
    int res1=run_test<3>("W(x,y,z)",tf1,coordf,"data_W.out");

    auto f2= W1d;
    void *tf2 = *(void**)(&f2);
    int res2=run_test<1>("W1d(x)",tf2,coordf,"./data_W1d.out");

    auto f3= Wp;
    void *tf3 = *(void**)(&f3);
    int res3=run_test<1>("W12(x)",tf3,coordf,"./data_W12.out");
    

    auto f4= W1;
    void *tf4 = *(void**)(&f4);
    int res4=run_test<1>("W1(x)",tf4,coordf,"./data_W1.out");

    auto f5= I_Wp;
    void *tf5 = *(void**)(&f5);
    int res5=run_test<2>("I_W12(a,b)",tf5,coordf,"./data_IW12.out");

    auto f6 = I_W1;
    void *tf6 = *(void**)(&f6);
    int res6=run_test<2>("I_W1(a,b)",tf6,coordf,"./data_IW1.out");

    return (res1 || res2 || res3 || res4 || res5 || res6);

}
