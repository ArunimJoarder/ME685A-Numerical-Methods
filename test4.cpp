#include"matrix.hpp"
#include"ode.hpp"
#include<cmath>

double f(double x, double y){
    return pow(x,2)*y - y;
}

// mat f_(double t, mat y){
//     mat f(2,1);
//     f.set(0, 0, -y[0][0]);
//     f.set(1, 0, y[1][0]);
//     return f;
// }

mat exact(double x_i, double x_f, int N){
    mat ans(1,N);
    double h = (x_f-x_i)/(N-1);
    for(int i = 0; i < N; i++){
        double x = x_i + h*i;
        ans.set(0,i, exp(pow(x,3)/3 - x));
    }
    return ans;
}

mat EulerI(double (*f)(double, double), double y_i, double x_i, double x_f, int N = 100){
    mat ans(1,N);
    double h = (x_f-x_i)/(N-1);

    double yi = y_i;
    double xi = x_i;
    for(int i = 0; i < N; i++){
        ans.set(0,i, yi);
        yi = yi/(1 - h*(pow(xi+h, 2) - 1));
        xi = xi + h;
    }
    return ans;
}

int main(){
    mat act = exact(0, 1, 5);
    act.disp();

    // mat rk4 = RK4(f, 1, 0, 1, 5);
    // rk4.disp();
        
    mat md = midPoint(f, 1, 0, 1, 5);
    md.disp();

    // mat pc = PC(f, 1, 0, 1, 5);
    // pc.disp();

    mat ei = EulerI(f, 1, 0, 1, 5);
    ei.disp();

    return 0;
}