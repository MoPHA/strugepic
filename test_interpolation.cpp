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
            std::cout << "Error is to large, Value is: " <<calc_val << " Should be:"<< ref_val << std::endl;  
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

    return 0;
}



int main(void){
    std::cout << "Running tests:" << std::endl;
    int retvals[2];
    std::string coordf[3]={"data_x.out","data_y.out","data_z.out"};
    retvals[0]=run_test("W1d",strugepic::W1d,coordf[0],"./data_W1d.out");
    retvals[1]=run_test("W",strugepic::W,coordf,"data_W.out");
    std::cout << retvals[0] << std::endl;
}
