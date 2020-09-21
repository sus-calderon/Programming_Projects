// <.h>
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
        //Variables:
        double enuc;    //Nuclear repulsion energy
        int norb;       //Number of atomic orbitals in file (water has 7, 1s for each H, and 1s, 2s, and 3 2p for O)
        Matrix S;       //The overlap integral is a 7x7 integral, and so s[i][j]=1.00000 if i=j=1
        Matrix T;
        Matrix V;
        Matrix core;    //My core Hamiltonian will be made up of the t and v (kinetic and potential)
        Matrix SOM;     //The Symmetric Orthogonal Overlap Matrix
        Vector TEI;     //Two electron repulsion integral
        Vector TEI_MO;  //Two electron repulsion integral in the MO basis
        Vector ioff;
        Matrix D;       //Density Matrix
        Matrix F;       //New Fock matrix formed from Density matrix and tei
        Matrix F_p;     //Orthogonalized Fock matrix
        Matrix C;       //New calculated coefficient matrix
        Matrix E_p;
        //Vector E_p;     //The Molecular Orbital energies
        int iter;
        int iter_max;
        double SCF;
        double tot_E;
        double old_SCF; //For my while loop, that tests for convergence
        Matrix old_D;
        //Variables for Smart MP2 algorithm
        Matrix tmp;
        Matrix X;
        Matrix Y;
        Matrix final_;



        //Functions
        void print_matrix(string mat_string, Matrix matrix);      
        void print_vector(string mat_string, Vector vect);      

        double read_enuc(const char *filename);
        int read_overlap(HartreeFock& hf, const char *filename);
        int read_kinetic(HartreeFock& hf, const char *filename);
        int read_potential(HartreeFock& hf, const char *filename);
        int build_core(HartreeFock& hf);
        int read_tei(HartreeFock& hf, const char *filename);

        int build_orthog(HartreeFock& hf);
        int build_density(HartreeFock& hf, int elec_num);
        
        int compute_SCF(HartreeFock& hf);
        int update_Fock(HartreeFock& hf);

        int transform_AO_2_MO(HartreeFock& hf); //This function is to transform from AO to MO basis
        int smart_transform_AO_2_MO(HartreeFock& hf); //N^5 smart algorithm
        int MP2_calc(HartreeFock& hf, int elec_num); //This fxn calculates MP2 energy
        int smart_MP2_calc(HartreeFock& hf, int elec_num); //This fxn calculates MP2 energy

        HartreeFock(const char *filename);
        ~HartreeFock();
};

