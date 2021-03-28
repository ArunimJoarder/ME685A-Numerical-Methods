#include<cmath>
#include"derivative.hpp"
#include"matrix.hpp"

const double pi = M_PI;
const double lambda = 9, alpha = 1, L = 1;
const vector<int> Ns = {101, 201};

// Values of del_t found for N= 101, and N = 201
    // del_t(N = 101) = 5.129000e-04, del_t(N = 201) = 5.178500e-04;


// Steady State distribution of T
double SS(double x){
    double a = pow(lambda,0.5);
    double T = sinh(a*(1-x))/sinh(a);
    return T;
}

// Analytical Solution of T as derived from asgn hints
double analyticalSol(double x, double t){
    double sum = 0;
    for(double n = 1.0; n <= 1000; n++){
        sum += 2*n*pi*sin(n*pi*x)/(n*n*pi*pi + lambda)*(1 - exp(-(n*n*pi*pi + lambda)*t));
    }
    return sum;
}

// Gradient of analytical solution as derived from asgn hints
double gradAnalyticalSol(double x, double t){
    double sum = 0;
    for(double n = 1; n <= 1000; n++){
        sum += 2*n*pi*n*pi*cos(n*pi*x)/(n*n*pi*pi + lambda)*(1 - exp(-(n*n*pi*pi + lambda)*t));
    }
    return sum;
}

// Returns coefficient matrix with FTCS approach
mat coeffMat(int N, double del_x, double del_t){
    mat a(N,N);
    
    // adding Dirichlet boundary conditions
    a.set(0,0, 1);
    a.set(N-1,N-1, 1);

    double Fo = alpha*del_t/(del_x*del_x);

    // setting coefficients for middle points
    for(int i = 1; i < N-1; i++){
        a.set(i,i-1, -Fo);
        a.set(i,i, (1+2*Fo+lambda*del_t));
        a.set(i,i+1, -Fo);
    }
    return a;
}

// main function to calculate variation of T with x and t
    // inputs initial conditions
vector<mat> FTCS_GS(mat T_i){
    
    // vector containing temperatures at all nodes for each time step
    vector<mat> temperatures;
    temperatures.push_back(T_i);

    int N = T_i.size()[0];
    
    // specify del_t as required
    double del_x = L/(N-1), del_t = 0.1;
    // if(N == 101)
    //     del_t = 5.129000e-04;
    // else if(N == 201)
    //     del_t = 5.178500e-04;


    // construct coefficient matrix
    mat A = coeffMat(N, del_x, del_t);
    // construct inverse using Gauss Seidel
    mat A_inv = inv(A);
    
    // calculate gradient of steady state solution at x = 0
    double ssGradient = d(SS, 0);
    cout << "\t2*steady_state_gradient = " << 2*ssGradient << endl;

    cout << "\n\tTemperature @ x = L/2\t\t\t\t\tTemperature Gradient @ x = 0\n\t  |  Numerical | Analytical |Steady-State|\t\t  |  Numerical  |   Analytical  | Steady-State|\n\t  |------------|------------|------------|\t\t  |-------------|---------------|-------------|\n";
    double gradient = INT16_MAX;
    int iter = 0;


    mat T = T_i, T_prev = T_i + 2;
    // check if gradient is greater than 2*steady_state_gradient
    while(norm(T-T_prev) > 1e-8 && iter < 500){
        // take b matrix as T at (n-1)th time step
        T_prev = temperatures.back();

        // calculate T at (n)th time step
        T = A_inv*T_prev;

        // add T at (n)th to the vector
        temperatures.push_back(T);

        // calculate gradient at x = 0
        gradient = (T[1][0] - T[0][0])/del_x;
        
        // num -> numerical answer, ana -> analytical answer
        cout << "\t  |  " << T[(N-1)/2][0] << "  |  " << analyticalSol(L/2, iter*del_t) << "  |  " << SS(L/2) << "  |" << "\t\t  |  " << gradient << "  |  " << gradAnalyticalSol(1e-10, iter*del_t) << "  |  " << ssGradient << "  |" << endl;
        
        // continue to next time step
        iter++;
    }

    cout << "\t  |------------|------------|------------|\t\t  |-------------|---------------|-------------|\n\n";
    // this statement was used find the required del_t
        // cout << "\t" << scientific << del_t*iter/20.0 << endl;
    return temperatures;
}

int main(){

    for(auto N : Ns){
        mat T_i(N,1);
        T_i.set(0,0,1);

        // do calculations for both values of N
        cout << fixed << "N = " << N << endl;
        vector<mat> ans = FTCS_GS(T_i);
    }

    return 0;
}