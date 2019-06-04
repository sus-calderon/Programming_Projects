//So this water.cc is my water molecule file.
//I will use my declared function that are declared in the .h file.
//I will name them with h20. which is a type of "Moleucle" and give them values.
//I will then have to compile my code to link the molecule.cc and molecule.h files 
// together into an executable program. (a .o file)
// .cc and such are raw source files we must compile as object file. We link object code to make them executable.
//
// Rewritten June 3, 2019
//

#include "molecule.h"

using namespace std;

int main(int argc, char *argv[])
{
    Molecule h20;

    h2o.natom = 3;
    h2o.charge = 0;
    h2o.zvals = new int[h2o.natom];
    h2o.geom = new double* [h2o.natom];
    for(int i=0; i < h2o.natom; i++)
        h2o.geom[i] = new double [3];

    h2o.zvals[0] = 8;
    h2o.geom[0][0] =  0.000000000000;
    h2o.geom[0][1] =  0.000000000000;
    h2o.geom[0][2] = -0.122368916506;
    h2o.zvals[1] = 1;
    h2o.geom[1][0] =  0.000000000000;
    h2o.geom[1][1] =  1.414995841403;
    h2o.geom[1][2] =  0.971041753535;
    h2o.zvals[2] = 1;
    h2o.geom[2][0] =  0.000000000000;
    h2o.geom[2][1] = -1.414995841403;
    h2o.geom[2][2] =  0.971041753535;

    h2o.print_geom();
    h2o.translate(5, 0, 0);
    h2o.print_geom();

    delete[] h2o.zvals;
    for(int i=0; i < h2o.natom; i++)
        delete[] h2o.geom[i];
    delete[] h2o.geom;

    return 0;
}
