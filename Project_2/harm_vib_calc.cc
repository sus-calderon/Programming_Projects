//This .cc file will use functions and variables I declared in a .h file 
//and definied in a corresponding .cc file
// Susana Calderon - June 17, 2019
// <vim-harb_calc.cc>
//

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <cmath>
#include "molecule.h"
#include "hessian.h"

using namespace std;

int main()
{
    //Read input file geometry
    Molecule mol("h2o_geom.txt", 0);
    cout << "Number of atoms: " << mol.natom << endl;
    cout << "Input Cartesian Coordinates: \n"; 
    mol.print_geom();
    cout << endl;

    //I now will read in the Cartesian Hessian Data (consists of 2nd derivatives of the E wrt atomic positions
    //Read in Hessian
    Hessian hess("h2o_hessian.txt"); 
    if (mol.natom != hess.natom) 
    {
        cout << "\n The hessian does not correspond to the geometry. (The number of atoms is not the same.)" << endl;
        exit(0);
    }
    //Read in and store the Hessian matrix values
    //Here I can use my new hessian.cc and hessian.h
    
    cout << "Input Hessian Matrix: \n";
    hess.print_H();
    cout << endl;

    return 0;
}

