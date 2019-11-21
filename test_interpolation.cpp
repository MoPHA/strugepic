#include<iostream>
#include"interpolation.hpp"
#include<iomanip>
#include <vector>
#include <fstream>
#include<cmath>
int main(void){

    std::cout << "Running tests:" << std::endl;
    std::vector<double> x_vals;
    std::vector<double> y_vals;
    std::vector<double> z_vals;
    std::vector<double> i_vals;
    std::vector<double> id_vals;

    std::ifstream xfile("data_x.out", std::ios::in);
    std::ifstream yfile("data_y.out", std::ios::in);
    std::ifstream zfile("data_z.out", std::ios::in);
    std::ifstream ifile("data_val.out", std::ios::in);
    std::ifstream idfile("data_vald.out",std::ios::in);
    
    double num = 0.0;
    
    if (!xfile.is_open()) {
        std::cerr << "There was a problem opening the input file!\n";
        exit(1);//exit or do additional error checking
    }
    if (!yfile.is_open()) {
        std::cerr << "There was a problem opening the input file!\n";
        exit(1);//exit or do additional error checking
    }
    if (!zfile.is_open()) {
        std::cerr << "There was a problem opening the input file!\n";
        exit(1);//exit or do additional error checking
    }
    if (!ifile.is_open()) {
        std::cerr << "There was a problem opening the input file!\n";
        exit(1);//exit or do additional error checking
    }
    if (!idfile.is_open()) {
        std::cerr << "There was a problem opening the input file!\n";
        exit(1);//exit or do additional error checking
    }


    while (xfile >> num) {
        x_vals.push_back(num);
    }
    while (yfile >> num) {
        y_vals.push_back(num);
    }
    while (zfile >> num) {

        z_vals.push_back(num);
    }
    while (ifile >> num) {
        i_vals.push_back(num);
    }
    while (idfile >> num) {
        id_vals.push_back(num);
    }
   
    int failed=0; 
    for(int i=0; i< (int) x_vals.size()-1;i++){
        if(i_vals[i] != 0){
            double relative_error = (strugepic::W(x_vals[i],y_vals[i],z_vals[i]) - i_vals[i] )/(i_vals[i]);
            if(fabs(relative_error) > 1.0e-10 && fabs(strugepic::W(x_vals[i],y_vals[i],z_vals[i]) - i_vals[i]) > 1.0e-9 ){
            std::cout << "error: value is " << strugepic::W(x_vals[i],y_vals[i],z_vals[i])<< std::setprecision (15)  
                << " should be: " << i_vals[i] << std::setprecision (15) << "  relative error is:" <<  relative_error <<std::endl;
            failed=1;
            }
        }
        else{
            if(strugepic::W(x_vals[i],y_vals[i],z_vals[i]) != 0){
            std::cout << "ERROR: value is " << strugepic::W(x_vals[i],y_vals[i],z_vals[i])<< std::setprecision (15)  
                << " Should be: " << i_vals[i] << std::setprecision (15) <<std::endl;
            std::cout << x_vals[i]  << " , " << y_vals[i] << " , " << z_vals[i] << "  " << i <<std::endl;
            failed=1;
            }
        }
    }
    for(int i=0; i< (int) x_vals.size()-1;i++){
        if(id_vals[i] != 0){
            double relative_error = (strugepic::W1d(x_vals[i]) - id_vals[i] )/(id_vals[i]);
            if(fabs(relative_error) > 1.0e-10 && fabs(strugepic::W1d(x_vals[i]) - id_vals[i]) > 1.0e-9 ){
            std::cout << "error: value is " << strugepic::W1d(x_vals[i])<< std::setprecision (15)  
                << " should be: " << id_vals[i] << std::setprecision (15) << "  relative error is:" <<  relative_error <<std::endl;
            failed=1;
        }
        }
    }
    if(failed==0){
        std::cout << "PASS" << std::endl;
    }
    return failed;

}
