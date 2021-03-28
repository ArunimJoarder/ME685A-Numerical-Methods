#include<cmath>
#include<iostream>
#include<vector>
#include"matrix.hpp"
#include"derivative.hpp"

#ifndef PICARD_H
#define PICARD_H
using namespace std;

double picard(double (*f)(double), double (*G)(double), double x0, double e = 0.01, int iter_max = 100){
    // main iterative loop
    int counter = 0;
    double x = x0;
    double x_prev = x + 2;
    double ea = e + 2;
    while((ea > e) && (fabs(d(G,x)) < 1) && (counter < iter_max)){
        counter++;
        x_prev = x;
        x = G(x);
        ea = fabs((x-x_prev)/x)*100;
    }
    if((fabs(d(G,x)) >= 1) && ea >= e){
        cout << "  Diverging Solution" << endl;
    }
    return x;
}
mat picard(mat (*f)(mat), mat (*G)(mat), mat x0, double e = 0.01, int iter_max = 100){
    // main iterative loop
    int counter = 0;
    mat x = x0;
    mat x_prev = x + 2;
    double ea = e + 2;
    while((ea > e) && det(d(G,x)) < 1 && (counter < iter_max)){
        counter++;
        x_prev = x;
        x = G(x);
        ea = norm(x-x_prev)/norm(x)*100;
    }
    if(det(d(f,x)) >= 1 && ea >= e){
        cout << "  Diverging Solution" << endl;
        return mat();
    }

    return x;
}
#endif