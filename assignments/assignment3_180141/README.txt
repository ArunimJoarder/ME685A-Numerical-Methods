The definitions of matrix functions are in the 'matrix.h' file.

The problem statement of the assignment is solved in the file 'asgn3.cpp'.

To compile program, type the following in your terminal;
    $ cd [path to where code is saved]
    $ g++ -std=c++11 asgn3.cpp -o asgn3     ## to set g++ compiler standard to C++11
    $ ./asgn3

OR, to run the executable file only;
    $ cd [path to where code is saved]
    $ ./asgn3

Explanation for change in Condition Number with change in Grid Fourier Number:
    We see that with increase in Grid Fourier Number the Condition Number also increases.
    This is because the the diagonal entry is -(1+2*Fo) and the sum of other row entires
    are 2*Fo, the difference in these (which describes Diagonal Dominance) is always 1.

    When the value of Fo is large, the difference is not much as comapred to the non-diagonal
    terms and thus the Diagonal Dominance is not to a high degree and is thus not easily invertible.

    
The code has been tested to be functioning on the mathserver.iitk.ac.in computer with above commands.