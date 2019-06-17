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
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
// A question I need to ask later is why did I have to include Eigen::Dynamic twice?
// Like why did it affect the order of my arrays and vectors, that was weird
//typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

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


    //Principal Moments of Inertia
    Matrix I(3,3);
    for(int i=0; i<mol.natom; i++)
    {
        mi = masses[(int) mol.zvals[i]];
        I(0,0) += mi*((mol.geom[i][1] * mol.geom[i][1]) + (mol.geom[i][2] * mol.geom[i][2]));
        I(1,1) += mi*((mol.geom[i][0] * mol.geom[i][0]) + (mol.geom[i][2] * mol.geom[i][2]));
        I(2,2) += mi*((mol.geom[i][0] * mol.geom[i][0]) + (mol.geom[i][1] * mol.geom[i][1]));
        I(0,1) += mi*(mol.geom[i][0])*mol.geom[i][1];
        I(0,2) += mi*(mol.geom[i][0])*mol.geom[i][2];
        I(1,2) += mi*(mol.geom[i][1])*mol.geom[i][2];
    }

    I(1,0) = I(0,1);
    I(2,0) = I(0,2);
    I(2,1) = I(1,2);

    cout << "\n Moment of Inertia Tensor (amu bohr^2): \n" << endl;
    cout << I << endl;

    Eigen::SelfAdjointEigenSolver<Matrix> solver(I);
    Matrix evecs = solver.eigenvectors();
    Matrix evals = solver.eigenvalues();

    cout << "\n Principal Moments of Inertia (amu * bohr^2); \n";
    cout << evals << endl;

    double conv = 0.529177249*0.529177249;
    cout << "\n Principal Moments of Inertia (amu * Angstrom^2): \n";
    cout << evals*conv << endl;

    double a=1.0e-16;
    double b=1.66054e-24;
    double conv2 = a*b;
    cout << "\n Principal Moments of Inertia (g * cm^2): \n";
    cout << evals*conv*conv2 << endl;

    

//    //For moment of Intertia, I want to find each element's momentf intertia tensor ( a 3x3 matrix)
//    Matrix I(3,3);  //Using the Eigen function I can make a matrix to store my inertia tensors in
//    for(int i=0; i<mol.natom; i++) 
//    {
//        mi = masses[(int) mol.zvals[i]];
//            I(0,0) += mi * (mol.geom[i][1]*mol.geom[i][1] + mol.geom[i][2]*mol.geom[i][2]);
//            I(1,1) += mi * (mol.geom[i][0]*mol.geom[i][0] + mol.geom[i][2]*mol.geom[i][2]);
//            I(2,2) += mi * (mol.geom[i][0]*mol.geom[i][0] + mol.geom[i][1]*mol.geom[i][1]);
//            I(0,1) += mi * mol.geom[i][0]*mol.geom[i][1];
//            I(0,2) += mi * mol.geom[i][0]*mol.geom[i][2];
//            I(1,2) += mi * mol.geom[i][1]*mol.geom[i][2];
//    }
    
//    I(1,0) = I(0,1);
//    I(2,0) = I(0,2);
//    I(2,1) = I(1,2);

//    cout << "\nMoment of inertia tensor (amu bohr^2):\n";
//    cout << I << endl;

//    // find the principal moments
//    Eigen::SelfAdjointEigenSolver<Matrix> solver(I);
//    Matrix evecs = solver.eigenvectors();
//    Matrix evals = solver.eigenvalues();

//    cout << "\nPrincipal moments of inertia (amu * bohr^2):\n";
//    cout << evals << endl;

//    double conv = 0.529177249 * 0.529177249;
//    cout << "\nPrincipal moments of inertia (amu * AA^2):\n";
//    cout << evals * conv << endl;

//    conv = 1.6605402E-24 * 0.529177249E-8 * 0.529177249E-8;
//    cout << "\nPrincipal moments of inertia (g * cm^2):\n";
//    cout << evals * conv << endl;

//    // classify the rotor
//    if(mol.natom == 2) cout << "\nMolecule is diatomic.\n";
//    else if(evals(0) < 1e-4) cout << "\nMolecule is linear.\n";
//    else if((fabs(evals(0) - evals(1)) < 1e-4) && (fabs(evals(1) - evals(2)) < 1e-4))
//        cout << "\nMolecule is a spherical top.\n";
//    else if((fabs(evals(0) - evals(1)) < 1e-4) && (fabs(evals(1) - evals(2)) > 1e-4))
//        cout << "\nMolecule is an oblate symmetric top.\n";
//    else if((fabs(evals(0) - evals(1)) > 1e-4) && (fabs(evals(1) - evals(2)) < 1e-4))
//        cout << "\nMolecule is a prolate symmetric top.\n";
//    else cout << "\nMolecule is an asymmetric top.\n";

    return 0;
}


