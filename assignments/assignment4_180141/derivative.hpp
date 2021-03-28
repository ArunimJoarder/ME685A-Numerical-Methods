#include<vector>
#include"matrix.hpp"

#ifndef DERIVATIVE_H
#define DERIVATIVE_H
using namespace std;

double d(double (*f)(double), double x){
    double e = 1e-3;
    double der = (f(x+e)-f(x-e))/(2*e);
    return der;
}

// double dd(double (*f)(double), double x){
//     double e = 1e-3;
//     double der = (f(x+e) - 2*f(x) + f(x-e))/(e*e);
//     // cout << 1 << " " << der << endl;
//     return der;
// }

double d(double (*f)(double, vector<double>), double x, vector<double> r){
    double e = 1e-3;
    double der = (f(x+e,r)-f(x-e,r))/(2*e);
    return der;
}

// double dd(double (*f)(double, vector<double>), double x, vector<double> r){
//     double e = 1e-3;
//     double der = (f(x+e,r) - 2*f(x,r) + f(x-e,r))/(e*e);
//     // cout << 1 << " " << der << endl;
//     return der;
// }

mat d(mat (*f)(mat), mat x){
    int n = f(x).size()[0];
    int m = x.size()[0];
    
    double E = 1e-3;
    mat e = E*eye(m);
    mat J(n,m);

    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            double g = ((f(x+e.col(j)) - f(x-e.col(j)))[i][0])/(2*E);
            J.set(i,j,g);
        } 
    }

    return J;
}

#endif