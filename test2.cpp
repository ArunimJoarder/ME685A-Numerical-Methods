#include<cmath>
#include<iostream>
#include<cstdio>
#include"integral.hpp"

using namespace std;

double f(double x){
    double a = exp(-x)*x;
    return a;
}
double ff(double x){
    double a = -(x+1)*exp(-x);
    return a;
}

int main(){
    int N[] = {5, 11, 41, 101};
    double a = 0.0, b = 3.0;
    double exact_sol = ff(b) - ff(a);
    printf("N\t   Exact      Trapezoidal Rule\t Simpson's Rule\t   Gauss 3-Point Quadrature\n");
    for(int i = 0; i < 4; i++){
        printf("%d\t%.8f\t %.8f\t   %.8f\t\t  %.8f\n", N[i], exact_sol, trapezoid(f,a,b,N[i]), simpson(f,a,b,N[i]), gauss3point(f,a,b,N[i]));
    }
    return 0;
}
