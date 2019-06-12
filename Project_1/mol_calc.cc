//This .cc file will use the functions and variables I declared in the .h file and defined in the .cc file.
//This will call back to what functions I need, and this will be used to actually analyze my molecules.
//Rewritten by Susana Calderon on June 10, 2019
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <string>
#include <cmath>
#include "molecule.h"
#include "masses.h"
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
    
    
    cout << "Number of atoms: " << mol.natom << endl;    //mol.natom is specific to mol file input
    
    
    cout << "Input Cartesian Coordinates:\n";    //Needs /n or endl so coord. aren't included in this line
    mol.print_geom();   //I can use my print_geom function to print out coordinates
    cout << endl;

    
    cout << "Interatomic Distances (bohr): \n";
    for(int i=0; i<mol.natom; i++)
        for(int j=0; j<i; j++)
            printf("%d %d %8.5f\n", i, j, mol.bond(i,j));
            //This prints out only unique values, where j<i
            //It solves for bond length then prints it, one by one, so I don't need to store data.
    cout << endl;

    
    cout << "Bond Angles: \n";
    for(int i=0; i<mol.natom; i++)
    {
        for(int j=0; j<i; j++)
        {
            for(int k=0; k<j; k++)
            {
                if(mol.bond(i,j)<4.0 && mol.bond(j,k)<4.0)  //don't need mol.bond(i,k) bc there's no real bond between them when j is the central atom
                {
                    printf("%2d- %2d- %2d %10.6f \n", i, j, k, mol.angle(i, j, k)*(180/acos(-1)) );
                }
            }
        }
    }
    cout << endl;


    cout << "Out-of-Plane Angles: \n";
    for(int i=0; i<mol.natom; i++) 
    {
        for(int k=0; k<mol.natom; k++) 
        {
            for(int j=0; j<mol.natom; j++) 
            {
                for(int l=0; l<j; l++) 
                {
                    if( i!=j && i!=k && i!=l && j!=k && j!=l && k!=l && mol.bond(i,k)<4.0 && mol.bond(k,j)<4.0 && mol.bond(k,l)<4.0 )
                        printf("%2d- %2d- %2d- %2d %10.6f \n", i, j, k, l, mol.oop(i,j,k,l)*(180.0/acos(-1.0)) );
                }
            }
        }
    }
    cout << endl;


    cout << "Torsional/Dihedral Angles: \n";
    for(int i=0; i<mol.natom; i++)
    {
        for(int j=0; j<i; j++)
        {
            for(int k=0; k<j; k++)
            {
                for(int l=0; l<k; l++)
                {
                    if( mol.bond(i,j)<4.0 && mol.bond(j,k)<4.0 && mol.bond(k,l)<4.0 )
                        printf("%2d- %2d- %2d- %2d %10.6f \n", i, j, k, l, mol.torsion(i,j,k,l)*(180.0/acos(-1.0)));
                }
            }
        }
    }
    cout << endl;


    //Find center of mass (COM); First initialize total mass variable
    double M = 0.0;
    for(int i=0; i<mol.natom; i++)  //We will use a loop to continually add the mass of each atom to total mass. (Use zval to ID atom)
    {
        M += masses[(int) mol.zvals[i]];    //This will go to masses.h and look for zvals position in array and use that value to add onto total masses.
    }

    double xcm = 0.0;
    double ycm = 0.0;
    double zcm = 0.0;
    double mi;

    for(int i=0; i<mol.natom; i++)
    {
        mi = masses[(int) mol.zvals[i]];
        xcm += mi*mol.geom[i][0];   //the formula asks to multiply mass by its x-position of that same atom and then add on for the next atom, it's mass and x position.
        ycm += mi*mol.geom[i][1];
        zcm += mi*mol.geom[i][2];
    }

    xcm /= M;
    ycm /= M;
    zcm /= M;

    printf("Molecular Center of Mass: %12.8f %12.8f %12.8f \n", xcm, ycm, zcm);
    cout << endl;
    mol.translate(-xcm, -ycm, -zcm);    //This translate my molecular geom to its COM as the origin (not the one given in the file)


    //For moment of Intertia, I want to find each element's momentf intertia tensor ( a 3x3 matrix)


    return 0;
}


