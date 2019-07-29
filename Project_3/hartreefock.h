// <hartreefock.h>
// My .h file will contain declaration of variables and functions.
// Written on July 15, 2019 by Susana Calderon
// Defined in hartreefock.cc

#include <string>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;  //Eigen::RowMajor makes it so that Eigen stores matrices by defaut in row-major order vs its actual default column-major order
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

using namespace std;

class HartreeFock     
{
    public:             //This mean that all the fxns and variables can be used by any code that makes use of an object called HartreeFock
        double enuc;    //Nuclear repulsion energy
        int norb;       //Number of atomic orbitals in file (water has 7, 1s for each H, and 1s, 2s, and 3 2p for O)
        Matrix S;     //The overlap integral is a 7x7 integral, and so s[i][j]=1.00000 if i=j=1
        Matrix T;
        Matrix V;
        Matrix core;  //My core Hamiltonian will be made up of the t and v (kinetic and potential)
        double *TEI;     //Two electron repulsion integral

        void print_matrix(string mat_string, Matrix matrix);      
        void read_enuc(const char *filename);
        Matrix read_overlap(Matrix oei_mat, const char *filename);
        Matrix read_kinetic(Matrix oei_mat, const char *filename);
        Matrix read_potential(Matrix oei_mat, const char *filename);
        Matrix build_core(Matrix t_mat, Matrix v_mat);
        void read_tei(double* tei_ary, const char *filename);
        void build_orthog(Matrix s_mat);
        void build_fock_guess(double** s_ortho, Matrix core_mat);

        HartreeFock(const char *filename);
        ~HartreeFock();
};

