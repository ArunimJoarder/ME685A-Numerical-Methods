#include"matrix.hpp"
#include<cmath>

// constructing coefficient matrix with DD BCs
mat DD(double Fo, int N){
    // initialise a (N x N) zero matrix
    mat A(N,N);

    // adding Dirichlet boundary conditions
    A.set(0,0, 1);
    A.set(N-1,N-1, 1);

    // setting coefficients for middle points
    for(int i = 1; i < N-1; i++){
        A.set(i,i-1, Fo);
        A.set(i,i, -(1+2*Fo));
        A.set(i,i+1, Fo);
    }
    return A;
}

// constructing coefficient matrix with DN BCs
mat DN(double Fo, int N){
    // initialise a (N x N) zero matrix
    mat A(N,N);

    // adding Dirichlet boundary conditions
    A.set(0,0, 1);

    // adding Newmann boundary conditions
    A.set(N-1,N-1, 1); A.set(N-1,N-2, -1);

    // setting coefficients for middle points
    for(int i = 1; i < N-1; i++){
        A.set(i,i-1, Fo);
        A.set(i,i, -(1+2*Fo));
        A.set(i,i+1, Fo);
    }
    return A;
}

// returns Condition Number of a matrix
double CN(mat A){
    // construct row-wise normalised matrix
    mat B(A);
    for(int i = 0; i < B.size()[0]; i++){
        double max_n = 0;
        for(int j = 0; j < B.size()[1]; j++){
            max_n = max(fabs(B[i][j]), max_n);
        }
        B.setRow(i, B.row(i)/max_n);
    }

    // find inverse of normalised matrix
    mat B_i = inv(B);

    printf("%.6f\t%.6f", row_sum_norm(B), row_sum_norm(B_i));
    // return product of norms of both matrices
    return row_sum_norm(B)*row_sum_norm(B_i);
}

int main(){
    auto GridNumbers = vector<double>{0.1, 0.25};
    int N = 5;

    printf("\tNOTE: B is the row normalised form of the coefficient matrix A\n\t\tThus, CN(A) = norm(B)*norm(B_inv)\n\n");
    printf("BC Type\t Fo\t  Norm(A)\tNorm(A_inv)\t Norm(B)       Norm(B_inv)\tCondition Number\n");
    for(auto Fo : GridNumbers){
        mat Add = DD(Fo, N);
        printf("  DD  \t%.2f\t %.6f\t %.6f\t", Fo, row_sum_norm(Add), row_sum_norm(inv(Add)));
        printf("\t   %.6f\n", CN(Add));
        mat Adn = DN(Fo, N);
        printf("  DN  \t%.2f\t %.6f\t %.6f\t", Fo, row_sum_norm(Adn), row_sum_norm(inv(Adn)));
        printf("\t   %.6f\n", CN(Adn));        
    }

    /*
        We see that with increase in Grid Fourier Number the Condition Number also increases.
        This is because the the diagonal entry is -(1+2*Fo) and the sum of other row entires
        are 2*Fo, the difference in these (which describes Diagonal Dominance) is always 1.

        When the value of Fo is large, the difference is not much as comapred to the non-diagonal
        terms and thus the Diagonal Dominance is not to a high degree and is thus not easily invertible.
    */
    return 0;
}