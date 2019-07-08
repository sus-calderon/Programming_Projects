//My .h file will contain declaration of variables and functions.
//I am writing this on June 17, 2019
//I will define them later in a .c file
// <hessian.h>

#include <string>

using namespace std;

class Hessian     //So this class is specific to reading in Hessians (which don't contain Z values 
{
    public:     //This mean that all the fxns and variables can be used by any code that makes use of an object called Hessian
        int natom;      //This should still exist since at the top of the Hessian values is the number of atoms
        int charge;
        double **2der;  //2der is a pointer-to-pointer-to-int array
        string point_group;  //Need to define a string varable although I don't know why it doesn't light up

        void print_2der();      //To print out Hessian matrix of 2nd derivatives

//May not this block for the Hessian
//        void rotate(double phi);
//        void translate(double x, double y, double z);
//        double bond(int atom1, int atom2);
//        double angle(int atom1, int atom2, int atom3);
//        double torsion(int atom1, int atom2, int atom3, int atom4);
//        double unit(int cart, int atom1, int atom2);
//        double oop(int atom1, int atom2, int atom3, int atom4);

        Hessian(const char *filename, int q);
        ~Hessian();
};

