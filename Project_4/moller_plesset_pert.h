// <moller_plesset_pert.h>
// My moller_plesset_pert.h file will contain declaration of variables and functions.
// Written on July 15, 2019 by Susana Calderon
// Defined in hartreefock.cc

#include <string>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"
#define INDEX(i,j) ((i>j) ? (((i)*((i)+1)/2)+(j)) :: (((j)*((j)+1)/2)+(i)))

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;  //Eigen::RowMajor makes it so that Eigen stores matrices by defaut in row-major order vs its actual default column-major order
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

using namespace std;

class MP2     
{
    public:
        //Variables:
        Vector TEI_AO;    //This exists in hf.cc as TEI
        Vector TEI_MO;

        //Functions
        //void print_matrix(string mat_string, Matrix matrix);      
        //void print_vector(string mat_string, Vector vect);      
        //I wonde if I can use the hf print functions

        int update_Fock(HartreeFock& hf);
        int transform_AO_2_MO(HartreeFock& hf, MP2& mp2);   //If I'm thinking about this correctly, it takes in two classes, my HF class so it can use the TEI array and MO coefficients and energies. Then it also has access to its own class for using TEI_AO and TEI_MO

        //void transform_2e_to_MO(HartreeFock& hf);



        MP2(const char *filename);
        ~MP2();
};

