/*
    ATTENTION
    This code is written in the C++11 standard, and may not compile for older compilers

    To run this code on the CCPC computer use the following command while in the code directory:
        $ g++ -std=c++11 asgn2.cpp -o asgn2
        $ ./asgn2
    
    This code has been tested on the mathsever.iitk.ac.in systems with the above commands
*/

#include<cmath>
#include<cstdio>
#include<iostream>

// Matrix datatype and associated function definitions are in 'matrix.h' file
#include"matrix.h"

using namespace std;

// Jacobian Definition
mat d(mat (*f)(mat), mat x){
    int n = f(x).size()[0];
    int m = x.size()[0];
    
    double E = 1e-9;
    mat e = E*eye(m);
    mat J(n,m);

    for(int j = 0; j < m; j++){
        mat g = (f(x+e.col(j)) - f(x-e.col(j)))/(2*E);
        J.setCol(j,g);
    } 

    return J;
}

// Newton-Raphson Iterations
mat newton(mat (*f)(mat),       // function for which roots are to be found
           mat x0,              // initial guess for root
           int iter_max = 100,  // default value of max number of iterations
           double e = 1e-5)     // default value for convergence criterion
{
    // main iterative loop
    int counter = 0;
    mat x = x0;
    mat x_prev = x + 2;
    double ea = e + 2;

    /* Definitions of det(X) and GE(A,b) are specified in the header file matrix.h */ 
    while((ea > e) && (fabs(det(d(f,x))) > 0.01) && (counter < iter_max)){
        // Increment number of iterations
        counter++;
        x_prev = x;
        // Update step
        x = x - GE(d(f,x), f(x));
        
        // Calculate relative error
        ea = norm(x-x_prev)/norm(x)*100;

        // Printing Statements for Table Entries
        printf("  %d\t[", counter);
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
        // In case soultion does not converge.
        cout << "  Diverging Solution" << endl;
        
        // returns 1x1 zero matrix
        return mat();
    }
    return x;
}

// Input function for root-finding
mat f(mat X){
    mat t(6,1);
    double Q1 = X[0][0], Q2 = X[1][0], Q3 = X[2][0], Q4 = X[3][0], Q5 = X[4][0], Q6 = X[5][0];
    double a, b, c, d, e, f;
    a = pow(Q1,2) + pow(Q2,2) + pow(Q3,2) - 14;         // F1
    b = 2*pow(Q2,2) + pow(Q3,2) + 2*pow(Q4,2) - 35;     // F2
    c = pow(Q3,2) + pow(Q4,2) - 2*pow(Q5,2) - 10;       // F3
    d = Q1 + Q5 - 3;                                    // F4
    e = Q4 + Q6 - 4;                                    // F5
    f = Q2 + Q6 - 3;                                    // F6
    
    // returns a column matrix with each ith row being Fi
    t.set(0,0, a);
    t.set(1,0, b);
    t.set(2,0, c);
    t.set(3,0, d);
    t.set(4,0, e);
    t.set(5,0, f);

    return t;
}

int main(){

    /*
        ATTENTION
        This code is written in the C++11 standard, and may not compile for older compilers

        To run this code on the CCPC computer use the following command while in the code directory:
            $ g++ -std=c++11 asgn2.cpp -o asgn2
            $ ./asgn2
        
        This code has been tested on the mathsever.iitk.ac.in systems with the above commands
    */

    // Table Headers
    printf("Index\t\t\tGuessed Values\t\t\t\t\tCorrections\t\t\t\t\tUpdates\n");
    
    // Initial Guess
    mat X0 = mat(arr2d {{5,5,5,5,5,5}}).T();
    
    mat x = newton(f,X0,9);

    // Print final soultion
    cout << "\nQ = [";
    for(int i = 0; i < x.size()[0]-1; i++){
        cout << x[i][0] << ", ";
    }
    cout << x[x.size()[0]-1][0] << "]\n";
    return 0;
    
    /*
        NOTE
        The solution converged before reaching 9 iterations, and therefore the code 
        doesn't generate table for the 9th iteration. This could be solved by removing
        convergence criteria for the newton-raphson method but I chose not to do that,
        as it meant the solution can misbehave and is not a good implementation of
        Newton-Raphson. Kindly take this into consideration.
    */
}