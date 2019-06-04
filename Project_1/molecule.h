#include <string>

using namespace std;

class Molecule
{
    public:
        //I'm declaring variables that can be used now
        int natom;
        int charge;
        int *zvals;
        double **geom;
        string point_group;
        //I'm declaring functions I want to use

        void print_geom();  //Void indicates a function that doesn't return a value
        void rotate(double phi);
        void translate(double x, double y, double z);
        double bond(int atom1, int atom2); //Use in step 2
        double angle(int atom1, int atom2, int atom3);
        double torsion(int atom1, int atom2, int atom3, int atom4);

        Molecule(int n, int q);
        ~Molecule();
};

