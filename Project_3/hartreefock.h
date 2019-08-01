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
        Matrix S;       //The overlap integral is a 7x7 integral, and so s[i][j]=1.00000 if i=j=1
        Matrix T;
        Matrix V;
        Matrix core;    //My core Hamiltonian will be made up of the t and v (kinetic and potential)
        Matrix SOM;     //The Symmetric Orthogonal Overlap Matrix
        Vector TEI;     //Two electron repulsion integral
        Matrix F_Guess; //Initial (guess) Fock Matrix
        Matrix MO_coef; //Inital Coefficient Matrix
        Matrix D;       //Density Matrix

        void print_matrix(string mat_string, Matrix matrix);      
        void print_vector(string mat_string, Vector vect);      

        void read_enuc(const char *filename);
        Matrix read_overlap(HartreeFock hf, const char *filename);
        Matrix read_kinetic(HartreeFock hf, const char *filename);
        Matrix read_potential(HartreeFock hf, const char *filename);
        Vector read_tei(HartreeFock hf, const char *filename);

        Matrix build_core(HartreeFock hf);
        Matrix build_orthog(HartreeFock hf);
        Matrix build_fock_guess(HartreeFock hf);
        Matrix build_MO_coef(HartreeFock hf);
        Matrix build_density(HartreeFock hf, int elec_num);

        HartreeFock(const char *filename);
        ~HartreeFock();
};

