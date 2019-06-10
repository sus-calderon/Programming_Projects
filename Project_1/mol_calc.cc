//This .cc file will use the functions and variables I declared in the .h file and defined in the .cc file.
//This will call back to what functions I need, and this will be used to actually analyze my molecules.
//Rewritten by Susana Calderon on June 20, 2019
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <string>
#include <cmath>
#include "molecule.h"
//#include "masses.h"
//#include "Eigen/Dense"
//#include "Eigen/Eigenvalues"
//#include "Eigen/Core"
//
//typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::RowMajor> Matrix;

using namespace std;

int main()
{
    //This will let my type in what file I want to be read in command line
    string filename;
    cout << "What molecule file do you want to read? \n";
    cin >> filename;

    //Now to call functions from molecule.cc
    Molecule mol(filename.c_str(), 0);  //I use filename..c_str() because the method to open files
                                        // takes in char[] when I have read in a string in the first
                                        // case. So i need to convert my string to characters.
                                        // "mol" is whatever file I have it read. It will be referred
                                        // to as mol in th rest of this code.
                                        // 0 = q for charge
    cout << "Number f atoms: " << mol.natom << endl;    //mol.natom is specific to mol file input
    cout << "Cartesian Coordinates: \n";    //Needs /n or endl so coord. aren't included in this line
    mol.print_geom();   //I can use my print_geom function to print out coordinates
    cout << endl;

    cout << "Interatomic Distances (bohr): \n";
    for(int i=0; i<mol.natom; i++)
        for(int j=0; j<i; j++)
            printf("%d %d %8.5f\n", i, j, mol.bond(i,j));
            //This prints out only unique values, where j<i
            //It solves for bond length then prints it, one by one, so I don't need to store data.
    cout << endl;

    return 0;
}


