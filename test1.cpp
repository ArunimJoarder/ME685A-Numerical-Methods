#include<cmath>
#include<iostream>
#include<vector>
#include"matrix.hpp"
#include"newton.hpp"
#include"false_pos.hpp"
#include"picard.hpp"

using namespace std;

// double f1(double x){
//     double a = pow(x,5) - 16.05*pow(x,4) + 88.75*pow(x,3) - 192.0375*pow(x,2) + 116.35*x + 31.6875;
//     return a;
// }

// double f2(double x, vector<double> r){
//     double a = pow(x,5) - 16.05*pow(x,4) + 88.75*pow(x,3) - 192.0375*pow(x,2) + 116.35*x + 31.6875;
//     for(auto xi : r){
//         a /= (x-xi);
//     }
//     return a;
// }

mat f3(mat X){
    mat t(6,1);
    double Q1 = X[0][0];
    double Q2 = X[1][0];
    double Q3 = X[2][0];
    double Q4 = X[3][0];
    double Q5 = X[4][0];
    double Q6 = X[5][0];
    double a = pow(Q1,2) + pow(Q2,2) + pow(Q3,2) - 14;
    double b = 2*pow(Q2,2) + pow(Q3,2) + 2*pow(Q4,2) - 35;
    double c = pow(Q3,2) + pow(Q4,2) - 2*pow(Q5,2) - 10;
    double d = Q1 + Q5 - 3;
    double e = Q4 + Q6 - 4;
    double f = Q2 + Q6 - 3;
    t.set(0,0, a);
    t.set(1,0, b);
    t.set(2,0, c);
    t.set(3,0, d);
    t.set(4,0, e);
    t.set(5,0, f);

    return t;
}

// mat f4(mat X){
//     mat t(2,1);
//     double x = X[0][0];
//     double y = X[1][0];
//     double a = f1(x);
//     double b = x+y;
//     t.set(0,0, a);
//     t.set(1,0, b);

//     return t;
// }

// double G1(double x){
//     double a = sqrt((f1(x) + 192.0375*pow(x,2))/192.0375);
//     return a;
// }

// mat G3(mat X){
//     mat t(2,1);
//     double x = X[0][0];
//     double y = X[1][0];
//     double a = sqrt((4-y)/2);
//     double b = sqrt((8-4*x));
//     t.set(0,0, a);
//     t.set(1,0, b);

//     return t;
// }

// mat G4(mat X){
//     mat t(2,1);
//     double x = X[0][0];
//     double y = X[1][0];
//     double a = sqrt((f1(x) + 192.0375*pow(x,2))/192.0375);
//     double b = -x;
//     t.set(0,0, a);
//     t.set(1,0, b);

//     return t;
// }

int main(){
    // // false position method
    // double xL = 0;
    // double xR = 6;
    // double xa = false_pos(f1,xL,xR);
    // cout << "(false_pos)\tx = " << xa << endl;

    // // bisection method
    // xa = bisection(f1,xL,xR);
    // cout << "(bisection)\tx = " << xa << endl;
    
    // // picard single root
    // double x0 = 0.0;
    // double xb = picard(f1, G1, x0);
    // cout << "(picard)\tx = " << xb << endl;

    // // newton single root
    // x0 = 0;
    // xb = newton(f1,x0);
    // cout << "(newton)\tx = " << xb << endl;

    // // newton multiple roots
    // x0 = 0;
    // vector<double> roots = newton(f2,x0);
    // cout << "(newton multiple roots)\troots = [";
    // for(int i = 0; i < roots.size()-1; i++){
    //     cout << roots[i] << ", ";
    // }
    // cout << roots.back() << "]\n";

    // // picard multivariable
    // mat X0 = mat(arr2d {{0, 0}}).T();
    // mat x = picard(f3, G3, X0);
    // cout << "(picard multivariate)\tx = [";
    // for(int i = 0; i < x.size()[0]-1; i++){
    //     cout << x[i][0] << ", ";
    // }
    // cout << x[x.size()[0]-1][0] << "]\n";

    // newton multivariable
    printf("Index\t\t\tGuessed Values\t\t\t\t\tCorrections\t\t\t\t\tUpdates\n");
    mat X0 = mat(arr2d {{5, 5, 5, 5, 5, 5}}).T();
    mat x = newton(f3,X0,1e-2,9);
    cout << "(newton multivariate)\tx = [";
    for(int i = 0; i < x.size()[0]-1; i++){
        cout << x[i][0] << ", ";
    }
    cout << x[x.size()[0]-1][0] << "]\n";
    return 0;
}