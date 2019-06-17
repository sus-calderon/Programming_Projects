//This .cc file will use functions and variables I declared in a .h file 
//and definied in a corresponding .cc file
// Susana Calderon - June 17, 2019
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
    string filename;
    cout << "What molecule file do you want to read? \n";
    cin >> filename;

    //We can use my molecule.h and molecule.cc files from the previous project so I can
    // call functions that read in and print the geometry. If I need any new functions
    // for this projects I can just add them in my .h and .cc file
    Molecule mol(filename.c_str(), 0);  // q = 0 for charge of my water
    cout << "Number of atoms: " << mo.natom << endl;
    cout << "Input Cartesian Coordinates: \n"; 
    mol.print_geom();
    cout endl;

    return 0;
}

