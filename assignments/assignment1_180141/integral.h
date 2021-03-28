#include<iostream>
#include<cmath>

#ifndef INTEGRAL_H
#define INTEGRAL_H
using namespace std;

double trapezoid(double (*f)(double), double a, double b, int N = 100){
    double del_x = (b-a)/(N-1);
    double sum = 0;
    for(int i = 2; i < N; i++){
        sum += f(a + (i-1)*del_x);
    }
    double integral = (f(a) + f(b) + 2*sum)*del_x/2;
    return integral;
}

double simpson(double (*f)(double), double a, double b, int N = 100){
    double del_x = (b-a)/(N-1);
    double sum_even = 0;
    for(int i = 2; i < N; i+=2){
        sum_even += f(a + (i-1)*del_x);
    }
    double sum_odd = 0;
    for(int i = 3; i < N; i+=2){
        sum_odd += f(a + (i-1)*del_x);
    }
    double integral = (f(a) + f(b) + 2*sum_odd + 4*sum_even)*del_x/3;
    return integral;
}

double gauss3point(double (*f)(double), double a, double b, int N = 10){
    double del_x = (b-a)/N;
    
    double integral = 0;
    for(int i = 0; i < N; i++){
        double w1 = 5/9.0, w2 = 8/9.0, w3 = 5/9.0;
        double z1 = -sqrt(3/5.0), z2 = 0, z3 = sqrt(3/5.0);
        double a1 = ((a + i*del_x) + (a + (i+1)*del_x))/2;
        double x1, x2, x3;
        x1 = del_x/2*z1 + a1;
        x2 = del_x/2*z2 + a1;
        x3 = del_x/2*z3 + a1;

        integral += (w1*f(x1) + w2*f(x2) + w3*f(x3))*del_x/2;
    }
    return integral;
}
#endif