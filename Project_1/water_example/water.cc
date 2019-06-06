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
    Molecule h2o("geom.dat",0);     //So this will read from the input file, what the zval is  

    h2o.print_geom();
    h2o.translate(5, 0, 0);
    h2o.print_geom();

    return 0;
}
