//So in .h files I declare variables and functions
//In .cc files I then define those variables and functions

#include "molecule.h"
#include <cstdio>

void Molecule::print_geom()
{
    for(int i=0; i<natom; i++)
        printf("%d %8.5f %8.5f %8.5f\n", zvals[i], geom[i][0], geom[i][1], geom[i][2]);
}

void Molecule::translate(double x, double y, double z)
{
    for(int i=0; i<natom; i++) {
        geom[i][0] += x;
        geom[i][1] += y;
        geom[i][2] += z;
    }
}

//The contructor and descructor do the work of allocating and deleting memory.
Molecule::Molecule(int n, int q)
{
    natom = n;
    charge = q;
    zvals = new int[natom];
    geom = new double* [natom];
    for(int i=0; i<natom; i++)
        geom[i] = new double[3];
}
Molecule::~Molecule()
{
    delete[] zvals;
    for(int i=0; i<natom; i++)
        delete[] geom[i];
    delete[] geom;
}
