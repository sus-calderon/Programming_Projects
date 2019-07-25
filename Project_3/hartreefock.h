// <hartreefock.h>
// My .h file will contain declaration of variables and functions.
// Written on July 15, 2019 by Susana Calderon
// Defined in hartreefock.cc

#include <string>

using namespace std;

class HartreeFock     
{
    public:             //This mean that all the fxns and variables can be used by any code that makes use of an object called HartreeFock
        double enuc;    //Nuclear repulsion energy
        int norb;       //Number of atomic orbitals in file (water has 7, 1s for each H, and 1s, 2s, and 3 2p for O)
        double **S;     //The overlap integral is a 7x7 integral, and so s[i][j]=1.00000 if i=j=1
        double **T;
        double **V;
        double **core;  //My core Hamiltonian will be made up of the t and v (kinetic and potential)
        double *R;     //Two electron repulsion integral

        void print_matrix(string mat_string, double** matrix);      
        void read_enuc(const char *filename);
        void read_oei(double** oei_mat, const char *filename);
        void build_core(double** t_mat, double** v_mat);
        void read_tei(double* tei_ary, const char *filename);
        void build_orthog(double** s_mat);

        HartreeFock(const char *filename);
        ~HartreeFock();
};

