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

using namespace std;

int main()
{
    //Read input file geometry
    //What's commented is for letting the user input what file to read
//    string filename;
//    cout << "What molecule file do you want to read? \n";
//    cin >> filename;
//    //We can use my molecule.h and molecule.cc files from the previous project so I can
//    // call functions that read in and print the geometry. If I need any new functions
//    // for this projects I can just add them in my .h and .cc file
//    Molecule mol(filename.c_str(), 0);  // q = 0 for charge of my water
    //These line will read a specific file I declare now in the code, and mol in the rest of the code
    //will relate back to atoms and coordinates to this specific file.
    Molecule mol("h2o_geom.txt", 0);
    cout << "Number of atoms: " << mol.natom << endl;
    cout << "Input Cartesian Coordinates: \n"; 
    mol.print_geom();
    cout << endl;

    //I now will read in the Cartesian Hessian Data (consists of 2nd derivatives of the E wrt atomic positions
    //Read in Hessian
    //Ok so this method is going to need changing. This reads input geometry of a molecule. 
    // But the hessian are derivatives, not coordinates persay and don't contain z-vals
    Molecule hess("h2o_hessian.txt", 0); 
    if (mol.natom != hess.natom) 
    {
        cout << "\n The hessian does not correspond to the geometry. (The number of atoms is not the same.)" << endl;
        exit(0);
    }
    //Read in and store the Hessian matrix values


    return 0;
}

