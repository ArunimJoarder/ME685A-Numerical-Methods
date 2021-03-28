#include<cmath>
#include<iostream>
#include<vector>
#include"matrix.hpp"
#include"derivative.hpp"

#ifndef NEWTON_H
#define NEWTON_H
using namespace std;

double newton(double (*f)(double), double x0, double e = 0.01, int iter_max = 100){
    // main iterative loop
    int counter = 0;
    double x = x0;
    double x_prev = x + 2;
    double ea = e + 2;
    while((ea > e) && (fabs(d(f,x)) > 1e-3) && (counter < iter_max)){
        counter++;
        x_prev = x;
        x = x - f(x)/d(f,x);
        ea = fabs((x-x_prev)/x)*100;
    }
    if(fabs(d(f,x)) >= 0.01 && ea >= e){
        cout << "  Diverging Solution" << endl;
    }
    return x;
}

vector<double> newton(double (*f)(double, vector<double>), double x0, double e = 0.01, int iter_max = 100){
    // parameters
    vector<double> r(0,0);
    bool flag = true;

    // main iterative loop
    while(flag){
        int counter = 0;
        double x = x0;
        double x_prev = x + 2;
        double ea = e + 2;
        while((ea > e) && (fabs(d(f,x,r)) > 1e-3) && (counter < iter_max)){
            counter++;
            x_prev = x;
            x = x - f(x,r)/d(f,x,r);
            ea = fabs((x-x_prev)/x)*100;
        }
        if(fabs(d(f,x,r)) >= 0.01 && ea >= e){
            cout << "  Diverging Solution" << endl;
            flag = false;
        }
        else if(counter == 0){
            flag = false;
        }
        else{
            r.push_back(x);
        }
    }
    return r;
}

mat newton(mat (*f)(mat), mat x0, double e = 0.01, int iter_max = 100){
    // main iterative loop
    int counter = 0;
    mat x = x0;
    mat x_prev = x + 2;
    double ea = e + 2;
    while((ea > e) && (fabs(det(d(f,x))) > 0.01) && (counter < iter_max)){
        counter++;
        x_prev = x;
        x = x - GE(d(f,x), f(x));
        ea = norm(x-x_prev)/norm(x)*100;

        // Printing Statements
        printf("%d\t[", counter);
        for(int i = 0; i < x_prev.size()[0]-1; i++){
            printf("%.4f,", x_prev[i][0]);
        }
        printf("%.4f]\t[", x_prev[x.size()[0]-1][0]);
        for(int i = 0; i < x.size()[0]-1; i++){
            printf("%.4f,", (x - x_prev)[i][0]);
        }
        printf("%.4f]\t[", (x - x_prev)[x.size()[0]-1][0]);
        for(int i = 0; i < x.size()[0]-1; i++){
            printf("%.4f,", x[i][0]);
        }
        printf("%.4f]\n", x[x.size()[0]-1][0]);
        //
    }
    if(fabs(det(d(f,x))) <= 0.01 && ea >= e){
        cout << "  Diverging Solution" << endl;
        return mat();
    }
    return x;
}
#endif