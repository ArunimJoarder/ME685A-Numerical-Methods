#include<iostream>
#include<vector>
#include<cmath>

#ifndef MATRIX_H
#define MATRIX_H
using namespace std;

typedef vector< vector<double> > arr2d;

typedef class Matrix{
    public:
        Matrix():r(1),c(1){
            vector<double> t(c,0.0);
            this->M.resize(r, t);
        }
        Matrix(const int m, const int n):r(m), c(n){
            vector<double> t(c,0.0);
            this->M.resize(r, t);
        }
        Matrix(arr2d A):r(A.size()), c(A[0].size()){
            vector<double> t(c,0.0);
            this->M.resize(r, t);
            for(int i = 0; i < r; i++){
                if(A[i].size() != A[0].size()){
                    cout << "Invalid matrix dimensions" << endl;
                    break;
                }
                for(int j = 0; j < c; j++){
                    this->M[i][j] = A[i][j];
                }
            }
        }
        Matrix(Matrix *A):r(A->size()[0]), c(A->size()[1]){
            M = A->getM();
        }

        // to set A[i,j] = val
        void set(int i, int j, double val){
            M[i][j] = val;
        }

        // to get 2D matrix in vector form
        arr2d getM(){
            return M;
        }

        // get (i,j)th value of matrix
        double get(int i, int j){
            return M[i][j];
        }

        vector<double> operator[](int i){
            return M[i];
        }
        void operator=(Matrix B){
            r = B.size()[0];
            c = B.size()[1];
            M = B.getM();
        }

        // transpose of matrix
        Matrix T(){
            Matrix t(c,r);
            for(int i = 0; i < r; i++){
                for(int j = 0; j < c; j++){
                    t.set(j,i, M[i][j]);
                }
            }
            return t;
        }

        // returns j-th column matrix 
        Matrix col(int j){
            Matrix t(r,1);
            if(j < c && j >= 0){
                for(int i = 0; i < r; i++){
                    t.set(i,0,this->M[i][j]);
                }
            }
            return t;
        }

        // returns i-th row matrix
        Matrix row(int i){
            Matrix t(1,c);
            if(i < r && i >= 0){
                for(int j = 0; j < c; j++){
                    t.set(0,j,this->M[i][j]);
                }
            }
            return t;
        }

        // sets value of i-th row of matrix
        void setRow(int i, Matrix row){
            if(row.size()[1] != c || row.size()[0] != 1){
                cout << "Row dimensions don't match" << endl;
                return;
            }
            for(int j = 0; j < c; j++){
                M[i][j] = row[0][j];
            }
        }

        // sets value of j-th column of matrix
        void setCol(int j, Matrix col){
            if(col.size()[0] != r || col.size()[1] != 1){
                cout << "Column dimensions don't match" << endl;
                return;
            }
            for(int i = 0; i < r; i++){
                M[i][j] = col[i][0];
            }
        }

        // displays values of matrix
        void disp(){
            // cout << endl;
            for(auto row :M){
                for(auto n : row){
                    printf("\t%.3f", n);
                }
                cout << endl;
            }
            cout << endl;
        }

        // returns the size of matrix as [no. of rows, no. of columns]
        vector<int> size(){
            vector<int> t(2,0);
            t[0] = r;
            t[1] = c;
            return t;
        }
    private:
        int r;
        int c;
        arr2d M;
}mat;

// returns identity matrix of size n
mat eye(int n){
    mat t(n,n);
    for(int i = 0; i < n; i++){
        t.set(i, i, 1);
    }
    return t;
}

// defines matrix addition
mat operator+(mat A, mat B){
    if(B.size() == A.size()){
        mat t(A);
        for(int i = 0; i < A.size()[0]; i++){
            for(int j = 0; j < A.size()[1]; j++){
                t.set(i,j,B[i][j]+A[i][j]);
            }
        }
        return t;
    }
    else{
        cout << "Incompatible matrix dimensions in addition" << endl;
        return A;
    }
}

// defines matrix and scalar addition
mat operator+(mat A, int k){
    mat t(A);
    for(int i = 0; i < A.size()[0]; i++){
        for(int j = 0; j < A.size()[1]; j++){
            t.set(i,j,k+A[i][j]);
        }
    }
    return t;
}
mat operator+(int k, mat A){
    return A+k;
}

// defines matrix multiplication
mat operator*(mat A, mat B){
    if(A.size()[1] == B.size()[0]){
        int m = A.size()[0];
        int n = A.size()[1];
        int p = B.size()[1];
        mat t(m,p);
        for(int i = 0; i < m; i++){
            for(int j = 0; j < p; j++){
                double sum = 0;
                for(int k = 0; k < n; k++){
                    sum += A[i][k]*B[k][j];
                }
                t.set(i,j,sum);
            }
        }
        return t;
    }
    else{
        cout << "Incompatible matrix dimensions in multiplication" << endl;
        return A;
    }
}

// defines matrix and scalar multiplication
mat operator*(double k, mat A){
    mat t(A);
    for(int i = 0; i < t.size()[0]; i++){
        for(int j = 0; j < t.size()[1]; j++){
            t.set(i,j,A[i][j]*k);
        }
    }
    return t;
}
mat operator*(mat A, double k){
    return k*A;
}

mat operator-(mat A){
    return -1*A;
}
mat operator-(mat A, mat B){
    return A+(-B);
}

// scalar divsion of matrices
mat operator/(mat A, double k){
    return A*(1.0/k);
}

// defines powers of matrices
mat operator^(mat A, double k){
    mat t = eye(A.size()[0]);
    for(int i = 0; i < k; i++){
        t = (t*A);
    }
    return t;
}

// norm for row and column matrices (magnitude of vectors)
double norm(mat A){
    if(A.size()[0] == 1){
        return sqrt((A*A.T())[0][0]);
    }
    else if(A.size()[1] == 1){
        return sqrt((A.T()*A)[0][0]);
    }
    else{
        // norm of matrices has not yet been included
        cout << "Not a vector (norm undefined)\n";
        return 0;
    }
}

// Check Diagonal dominance of matrix
bool checkDD(mat A){
    if(A.size()[0] != A.size()[1]){
        cout << "Not Diagonally Dominant (Not Square)\n";
        return false;
    }
    int n = A.size()[0];
    int check = 0;
    for(int i = 0; i < n; i++){
        int aii = 0;
        int aij = 0;
        for(int j = 0; j < n; j++){
            if(i == j){
                aii = fabs(A[i][j]);
            }
            else{
                aij += fabs(A[i][j]);
            }
        }
        if(aii == aij) check++;
        else if(aii > aij) check += 2;
        else{
            cout << "Not Diagonally Dominant\n";
            return false;
        }
    }
    if(check > n)
        return true;
    cout << "Not Diagonally Dominant\n";
    return false;
}

// Gauss Seidel Method of Solving Linear Equations
mat GS(mat A, mat b, double e1 = 0.01, double e2 = 0.001, double al = 1.5, int iter_max = 100){
    /* A and b are taken such that
        x = inv(A)*b
    */
    if(!checkDD(A))
        return mat();
    if(A.size()[0] != b.size()[0]){
        cout << "Dimensions of matrices don't match\n";
        return mat();
    }
    
    // parameters
    int n = A.size()[0];
    
    // initial guess
    mat x0 = b;
    for(int i = 0; i < b.size()[0]; i++){
        for(int j = 0; j < b.size()[1]; j++){
            x0.set(i,j,b[i][j]/A[i][i]);
        }
    }

    mat x = x0;
    mat x_prev = x;
    mat x_new = x;
    double ea1 = e1 + 2;
    double ea2 = e2 + 2;
    double counter = 0;

    while(ea1 > e1 || ea2 > e2 && counter < iter_max){
        counter++;
        x_prev = x;
        for(int i = 0; i < n; i++){
            for(int k = 0; k < b.size()[1]; k++){
                double sum = 0;
                for(int j = 0; j < n; j++){
                    sum += A[i][j]*x[j][k];
                }
                double a = (b[i][k] - sum + A[i][i]*x[i][k])/A[i][i];
                x_new.set(i,k,a);
                x.set(i,k,al*x_new[i][k] + (1-al)*x_prev[i][k]);
            }
        }
        ea1 = norm(x-x_prev)/sqrt(n);
        ea2 = norm(A*x-b);
    }
    return x;
}

// Gauss Elimination Method of solving Linear Equations
mat GE(mat A, mat b){
    /* A and b are taken such that
        x = inv(A)*b
    */
    if(A.size()[0] != A.size()[1]){
        cout << "Not Square\n";
        return mat();
    }
    if(A.size()[0] != b.size()[0]){
        cout << "Dimensions of matrices don't match\n";
        return mat();
    }

    // parameters
    int n = A.size()[0];
    
    // forward elimination
    for(int i = 0; i < n; i++){
        for(int j = i+1; j < n; j++){
            if(A[i][i] == 0 && A[j][i] != 0){
                mat r = A.row(i);
                A.setRow(i, A.row(j));
                A.setRow(j,r);
                mat br = b.row(i);
                b.setRow(i, b.row(j));
                b.setRow(j,br);
            }
            double coef = A[j][i]/A[i][i];
            
            mat r = A.row(j) - coef*A.row(i);
            A.setRow(j, r);

            mat br = b.row(j) - coef*b.row(i);
            b.setRow(j, br);
        }
    }

    // back substitution
    mat x(b);
    for(int i = n-1; i >= 0; i--){
        for(int k = 0; k < b.size()[1]; k++){
            double sum = 0;
            for(int j = i+1; j < n; j++){
                sum += A[i][j]*x[j][k];
            }
            double temp = (b[i][k] - sum)/A[i][i];
            x.set(i,k,temp);
        }
    }
    return x;
}


// Determinant of Matrix through Gauss Elimination method
double det(mat A){
    if(A.size()[0] != A.size()[1]){
        cout << "Not Square\n";
        return 0;
    }

    // parameters
    int n = A.size()[0];
    
    // forward elimination
    for(int i = 0; i < n; i++){
        for(int j = i+1; j < n; j++){
            if(A[i][i] == 0 && A[j][i] != 0){
                mat r = A.row(i);
                A.setRow(i, A.row(j));
                A.setRow(j,r);
            }
            double coef = A[j][i]/A[i][i];
            
            mat r = A.row(j) - coef*A.row(i);
            A.setRow(j, r);
        }
    }
    double prod = 1;
    for(int i = 0; i < n; i++){
        prod *= A[i][i];
    }
    return prod;
}

// Inverse of matrix through Gauss Jordan Method
mat inv(mat A){
    if(A.size()[0] != A.size()[1]){
        cout << "Not Square\n";
        return mat();
    }

    // parameters
    int n = A.size()[0];
    mat A_inv = eye(n);
    
    // forward elimination
    for(int i = 0; i < n; i++){
        for(int j = i+1; j < n; j++){
            if(A[i][i] == 0 && A[j][i] != 0){
                mat r = A.row(i);
                A.setRow(i, A.row(j));
                A.setRow(j,r);
                mat ir = A_inv.row(i);
                A_inv.setRow(i, A_inv.row(j));
                A_inv.setRow(j,ir);
            }
            double coef = A[j][i]/A[i][i];
            
            mat r = A.row(j) - coef*A.row(i);
            A.setRow(j, r);

            mat ir = A_inv.row(j) - coef*A_inv.row(i);
            A_inv.setRow(j, ir);
        }
    }
    // backward elimination
    for(int i = n-1; i >= 0; i--){
        for(int j = i-1; j >= 0; j--){
            if(A[i][i] == 0 && A[j][i] != 0){
                cout << "Division by 0\n";
                return mat();
            }
            double coef = A[j][i]/A[i][i];
            
            mat r = A.row(j) - coef*A.row(i);
            A.setRow(j, r);

            mat ir = A_inv.row(j) - coef*A_inv.row(i);
            A_inv.setRow(j, ir);
        }
    }
    // normalising diagonal terms
    for(int i = n-1; i >= 0; i--){
        if(A[i][i] == 0){
            cout << "Division by 0\n";
            return 0;
        }
        double coef = A[i][i];
        
        A.set(i,i, 1);

        mat ir = A_inv.row(i)/coef;
        A_inv.setRow(i, ir);
    }
    return A_inv;
}

// returns row-sum norm of a matrix
double row_sum_norm(mat A){
    double max_sum = 0;
    for(auto row : A.getM()){
        double sum = 0;
        for(auto aij : row){
            sum += fabs(aij);
        }
        // finds max among summmation of absolute entries of all rows
        max_sum = max(max_sum,sum);
    }
    return max_sum;
}

#endif