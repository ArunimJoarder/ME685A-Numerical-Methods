The definitions of matrix functions are in the 'matrix.hpp' file.
The definitions of derivative functions are in the 'derivative.hpp' file.


The problem statement of the assignment is solved in the file 'asgn4.cpp'.

To compile program, type the following in your terminal;
    $ cd [path to where code is saved]
    $ g++ -std=c++11 asgn4.cpp -o asgn4     ## to set g++ compiler standard to C++11
    $ ./asgn4

OR, to run the executable file only;
    $ cd [path to where code is saved]
    $ ./asgn4

NOTE:
The given analytical solution for T does not satisfy boundary conditions of T(0,t) = 1
In the solution given in the hints, the value of T(0,t) = 0 for all t.
This leads to a disparity in the solutions derived by analytical method, and the solution found by the numerical method.

The code has been tested to be functioning on the mathserver.iitk.ac.in computer with above commands.