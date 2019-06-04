//So in .h files I declare variables and functions
//In .cc files I then define those variables and functions

#include "molecule.h"
#include <cstdio>       //this is the C Standard Iput and Output library
//So in order to read from input files I need to include specific libraries that can do such
#include <iostream>     //header that defines the standard input/output stream objects
#include <fstream>      //input output stream to operate on files
#include <iomanip>      //header that provides parametric manipulators
#include <cassert>      //provides assert and standard debugging tool

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
Molecule::Molecule(const char *filename, int q)
{
    charge = q;

    //Open file
    std::ifstream is(filename);
    assert(is.good());
    //Read the number of atoms in the first line
    is >> natom;
    //Allocate space
    zvals = new int[natom];
    geom = new double* [natom];
    for(int i=0; i<natom; i++)
        geom[i] = new double[3];

    for(unsigned int i=0; i<natom; i++)
        is >> zvals[i] >> geom[i][0] >> geom[i][1] >> geom[i][2];   //specify the order in which the file is read into the input stream

    is.close();
}

Molecule::~Molecule()
{
    delete[] zvals;
    for(int i=0; i<natom; i++)
        delete[] geom[i];
    delete[] geom;
}
