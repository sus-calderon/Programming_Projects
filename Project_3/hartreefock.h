// <hartreefock.h>
// My .h file will contain declaration of variables and functions.
// Written on July 15, 2019 by Susana Calderon
// Defined in hartreefock.cc

#include <string>

using namespace std;

class HartreeFock     
{
    public:             //This mean that all the fxns and variables can be used by any code that makes use of an object called HartreeFock
        int norb;       //Number of atomic orbitals in file (water has 7, 1s for each H, and 1s, 2s, and 3 2p for O)
        double **S;     //The overlap integral is a 7x7 integral, and so s[i][j]=1.00000 if i=j=1
                        //My program is going to have to read in those values as i and j for my integral
        double **T;
        double **V;
        double **core;  //My core Hamiltonian will be made up of the t and v (kinetic and potential)

        void print_matrix(string mat_string, double** matrix);      
        void read_oei(double** oei_mat, const char *filename);
        void build_core(double** t_mat, double** v_mat);

        HartreeFock(const char *filename);
        ~HartreeFock();
};

