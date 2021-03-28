#include<iostream>
#include<cmath>
#include"matrix.hpp"

#ifndef ODE_H
#define ODE_H
using namespace std;

mat EulerE(double (*f)(double, double), double y_i, double x_i, double x_f, int N = 101){
    mat ans(1,N);
    double h = (x_f-x_i)/(N-1);

    double yi = y_i;
    double xi = x_i;
    for(int i = 0; i < N; i++){
        ans.set(0,i, yi);
        yi = yi + h*f(xi,yi);
        xi = xi + h;
    }
    return ans;
}

mat midPoint(double (*f)(double, double), double y_i, double x_i, double x_f, int N = 101){
    mat ans(1,N);
    double h = (x_f-x_i)/(N-1);

    double yi = y_i;
    double xi = x_i;
    for(int i = 0; i < N; i++){
        ans.set(0,i, yi);
        double yi_2 = yi + h/2*f(xi,yi);
        yi = yi + h*f(xi+h/2, yi_2);
        xi = xi + h;
    }
    return ans;
}

mat PC(double (*f)(double, double), double y_i, double x_i, double x_f, int N = 101){
    mat ans(1,N);
    double h = (x_f-x_i)/(N-1);

    double yi = y_i;
    double xi = x_i;
    for(int i = 0; i < N; i++){
        ans.set(0,i, yi);
        double yi_ = yi + h*f(xi,yi);
        yi = yi + h/2*(f(xi+h, yi_) + f(xi,yi));
        xi = xi + h;
    }
    return ans;
}

mat RK4(double (*f)(double, double), double y_i, double x_i, double x_f, int N = 101){
    mat ans(1,N);
    double h = (x_f-x_i)/(N-1);

    double yi = y_i;
    double xi = x_i;
    for(int i = 0; i < N; i++){
        ans.set(0,i, yi);
        xi = xi + h;
        double K1 = f(xi, yi);
        double K2 = f(xi + h/2, yi + K1*h/2);
        double K3 = f(xi + h/2, yi + K2*h/2);
        double K4 = f(xi + h, yi + K3*h);
        yi = yi + h/6*(K1 + 2*K2 + 2*K3 + K4);
    }
    return ans;
}

mat RK4(mat (*f)(double, mat), mat y_i, double t_i, double t_f, int N = 101){
    int n = f(t_i, y_i).size()[0];
    mat ans(n,N);
    double h = (t_f-t_i)/(N-1);

    mat yi = y_i;
    double ti = t_i;
    for(int i = 0; i < N; i++){
        ans.setCol(i, yi);
        ti = ti + h;
        mat K1 = f(ti, yi);
        mat K2 = f(ti + h/2, yi + K1*h/2);
        mat K3 = f(ti + h/2, yi + K2*h/2);
        mat K4 = f(ti + h, yi + K3*h);
        yi = yi + h/6*(K1 + 2*K2 + 2*K3 + K4);
    }
    return ans;
}
#endif