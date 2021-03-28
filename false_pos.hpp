#include<cmath>
#include<iostream>

#ifndef FALSE_POS_H
#define FALSE_POS_H
using namespace std;

double false_pos(double (*f)(double), double XL, double XR, double e = 0.01, int iter_max = 100){
    // False Position Method
    int counter = 0;
    double xL = XL;
    double xR = XR;
    double xr = xL;
    double x_prev = xr + 2;
    double ea = e + 2;
    while (fabs(ea) > e && counter < iter_max){
        counter++;
        if(f(xL)*f(xR) < 0){
            x_prev = xr;
            xr = (f(xR)*xL - f(xL)*xR)/(f(xR) - f(xL));
            if(f(xr)*f(xL) < 0)
                xR = xr;
            else if(f(xr)*f(xR) < 0)
                xL = xr;
        }
        ea = (xr - x_prev)/xr*100;
    }
    return xr;
}

double bisection(double (*f)(double), double XL, double XR, double e = 0.01, int iter_max = 100){
    // Bisection Method
    int counter = 0;
    double xL = XL;
    double xR = XR;
    double xr = xL;
    double x_prev = xr + 2;
    double ea = e + 2;
    while (fabs(ea) > e && counter < iter_max){
        counter++;
        if(f(xL)*f(xR) < 0){
            x_prev = xr;
            xr = (xL + xR)/2;
            if(f(xr)*f(xL) < 0)
                xR = xr;
            else if(f(xr)*f(xR) < 0)
                xL = xr;
        }
        ea = (xr - x_prev)/xr*100;
    }
    return xr;
}

#endif