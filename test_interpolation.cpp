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

// 0 -> test Ok
// 1 -> Numeric test failed,
// 2 -> could not open input
int run_test(std::string test_name,double (*test_func)(double),std::string arg_data_file,std::string ref_data_file  ){
    
    double num=0.0;
    std::vector<double> coord_vals;
    std::vector<double> ref_vals;
    std::ifstream coordfile(arg_data_file, std::ios::in);
    std::ifstream reffile(ref_data_file, std::ios::in);
    if (!coordfile.is_open() || !reffile.is_open()) {
        return 2;
    }
    
    while (coordfile >> num) {
        coord_vals.push_back(num);
    }
    while (reffile >> num) {
        ref_vals.push_back(num);
    }
    int cr;
    for(int i=0; i< (int) coord_vals.size()-1;i++){
     
        cr=compare_res(ref_vals[i],test_func(coord_vals[i]));
        if(cr !=0){
            return 1;
        }
    }

    return 0;
}; 
int run_test(std::string test_name,double (*test_func)(double,double,double),std::string *arg_data_files,std::string ref_data_file  ){

    double num=0.0;
    std::vector<double> coord_vals[3];
    std::vector<double> ref_vals;
    std::ifstream coordfiles[3];
    coordfiles[0].open(arg_data_files[0], std::ios::in);
    coordfiles[1].open(arg_data_files[1], std::ios::in);
    coordfiles[2].open(arg_data_files[2], std::ios::in);
    std::ifstream reffile(ref_data_file, std::ios::in);
    
    if (!coordfiles[0].is_open() || !coordfiles[1].is_open() || !coordfiles[2].is_open() || !reffile.is_open() ) {
        return 2;
    }

    for(int i=0; i<3; i++){
        while (coordfiles[i] >> num) {
            coord_vals[i].push_back(num);
        }
    } 

        while (reffile >> num) {
            ref_vals.push_back(num);
        }

    int cr;
    for(int i=0; i< (int) ref_vals.size()-1;i++){
     
        cr=compare_res(ref_vals[i],test_func(coord_vals[0][i],coord_vals[1][i],coord_vals[2][i]));
        if(cr !=0){
            return 1;
        }
    }
    
    return 0;
}



int main(void){
    int retval=0;
    std::cout << "Running tests:" << std::endl;

    std::string coordf[3]={"data_x.out","data_y.out","data_z.out"};
    std::cout << " W(x,y,z):";
    std::cerr << " W(x,y,z):" << std::endl;
    retval=run_test("W",strugepic::W,coordf,"data_W.out");
    retval? std::cout << " FAIL" : std::cout << " PASS";
    std::cout << std::endl;

    std::cout << " W1d(x):" ;
    std::cerr << " W1d(x):" << std::endl;
    retval=run_test("W1d",strugepic::W1d,coordf[0],"./data_W1d.out");
    retval? std::cout << " FAIL" : std::cout << " PASS";
    std::cout << std::endl;
}
