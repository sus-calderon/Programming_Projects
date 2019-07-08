//My .h file will contain declaration of variables and functions.
//I am rewriting this on June 8, 2019
//I will define them later in a .c file

#include <string>

using namespace std;

class Molecule      //So the name of our class will be Molecule
{
    public:     //This mean that all the fxns and variables can be used by any code
                // that makes use of an object called Molecule
        int natom;
        int charge;
        double *zvals;     //zvals will be an array so it needs a pointer * (pointer to int)
        double **geom;  //geom is a pointer-to-pointer-to-int array
        string point_group;  //Need to define a string varable although I don't know why it doesn't light up

        void print_geom();
        void rotate(double phi);
        void translate(double x, double y, double z);
        double bond(int atom1, int atom2);
        double angle(int atom1, int atom2, int atom3);
        double torsion(int atom1, int atom2, int atom3, int atom4);
        double unit(int cart, int atom1, int atom2);
        double oop(int atom1, int atom2, int atom3, int atom4);

        Molecule(const char *filename, int q);
        ~Molecule();
};

