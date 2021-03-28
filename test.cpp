#include"matrix.hpp"
#include<cmath>

int main(){
    mat A = arr2d{{1, 0, 0, 0},
                  {1, -2, 1, 0},
                  {0, 1, -2, 1},
                  {0, 0, -1, 1}};

    mat L = arr2d{{1, 0, 0, 0},
                  {1, 1, 0, 0},
                  {0, -0.5, 1, 0},
                  {0, 0, 2/3.0, 1}};

    mat U_ = arr2d{{1, 0, 0, 0},
                   {0, -2, 0, 0},
                   {0, 0, -1.5, 0},
                   {0, 0, 0, 1/3.0}};

    // A.disp();
    // mat b = arr2d{{1},
    //               {2},
    //               {0}};
    // cout << CN(A) << endl;

    mat B = inv(L*U_)*A;
    // B.disp();

    cout << CN(B) << endl;
    // mat x = GE(A,b);
    // x.disp();
    // mat A_ = GramSchimdt(A,b).first;
    // mat b_ = GramSchimdt(A,b).second;
    // A_.disp();
    // b_.disp();
    // cout << CN(A_) << endl;
    // x = GE(A_,b_);
    // x.disp();
    return 0;
}