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
#include "masses.h"

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

using namespace std;

int main()
{
    //Read input file geometry
    Molecule mol("h2o_geom.txt", 0);
    cout << "Number of atoms: " << mol.natom << endl;
    cout << endl;
    cout << "Input Cartesian Coordinates: \n"; 
    mol.print_geom();
    cout << endl;

    //Read in Cartesian Hessian Data (consists of 2nd derivatives of the E wrt atomic positions)
    //Here I can use my new hessian.cc and hessian.h
    Hessian hess("h2o_hessian.txt"); 
    if (mol.natom != hess.natom) 
    {
        cout << "\n The hessian does not correspond to the geometry. (The number of atoms is not the same.)" << endl;
        exit(0);
    }

    //Read in and store the Hessian matrix values
    cout << "Input Hessian Matrix: \n";
    hess.print_H();
    cout << endl;

    //Mass Weighted Hessian Matrix
    //I need to find the mass of my atoms using molecule.h and masses.h
    double m[mol.natom];
    for(int i=0; i<mol.natom; i++)
        m[i] = masses[(int) mol.zvals[i]];

    cout << "Weight of atoms (amu): \n";
    for(int i=0; i<mol.natom; i++)
        printf("%14.6f \n", m[i]);
    cout << endl;
    
    for(int i=0; i<3*mol.natom; i++) {
        for(int j=0; j<3*mol.natom; j++) {
            hess.H[i][j]/=sqrt(m[i/3]*m[j/3]);
        }
    }
    cout << "Mass Weighted Hessian Matrix: \n";
    hess.print_H();
    cout << endl;

    //Diagonalize the Hessian using the Eigen package
    Matrix F(3*mol.natom, 3*mol.natom);
    for(int i=0; i<3*mol.natom; i++) {
        for(int j=0; j<3*mol.natom; j++) {
            F(i,j) = hess.H[i][j];
        }
    }
    //cout << F << endl;
    //cout << endl;

    //This will solve for the eigenvalues
    Eigen::SelfAdjointEigenSolver<Matrix> solver(F);
    Matrix evecs = solver.eigenvectors();
    Matrix evals = solver.eigenvalues();

    cout << "Eigenvalues of the Hessian Matrix (hartee/amu*bohr^2): \n";
    cout << evals << endl;
    
    cout << endl;
    int size = evals.size();
    double w[size];
    
    //Eigenvalues are in hartee/amu*bohr^2
    //Conversion values obtained from NIST database
    //First conversion is from Hartee to Joules who's SI units are kg*m^2*s^-2
    double conv = 1.0;
    conv /= 2.293712278e17; // 1 Joule is equal to 2.293..x10^17 Hartrees

    //Next conversion is amu to kg
    conv /= 1.660539066e-27;  // 1amu = 1.6611x10^-27kg

    //Next conversion is bohr to m
    conv /= (5.291772109e-11)*(5.291772109e-11);  //This is bohr^2 & 1bohr = 5.29..x10^-11m

    double c = 2.99792458e10;    // Speed of light in cm/s
    double pi = 3.141592653;    //Value of pi

    for(int i=0; i<size; i++) 
        w[i] = (sqrt(conv*evals(i)))/(2*pi*c);
        //I also need to divide by 2*pi because when finding frequency, I want to know
        // how many cycles per second I have a full cycle(circle) is 2*pi
    
    cout << "Harmonic Vibrational Frequencies (cm^-1): \n";
    for(int i=0; i<size; i++)
        printf("%12.6f \n", w[i]);
    cout << endl;

    return 0;
}

