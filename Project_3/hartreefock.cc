// <hartreefock.cc>
// Written on July 16, 2019 by Susana Calderon
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cassert>
#include <cmath>
#include "hartreefock.h"
#include "molecule.h"
#include "electrons.h"

//Include Eigen package for easy diagonalization
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;  //Eigen::RowMajor makes it so that Eigen stores matrices by defaut in row-major order vs its actual default column-major order
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;


//Allocate memory for matrices after reading in norb = # of AOs 
HartreeFock::HartreeFock(const char *filename)
{
    //Open File here
    std::ifstream is(filename);
    assert(is.good());
   
    //This snippet of codes goes into the last line of file, sets the cursor to the last line,
    //and makes sure there is something to read in that last line.
    is.seekg(0,std::ios_base::end);         //Start at last line
    char ch = ' ';                          //Init ch not equal to '\n'
    while(ch != '\n'){
        is.seekg(-2,std::ios_base::cur);    //2 steps back, meaning it doesn't check last character
        if((int)is.tellg() <= 0 ){          //If it passed the start of the file, this is the start of the line
            is.seekg(0);
            break;
        }
        is.get(ch); //Checks the next character
    }

    //Grab first value of last line and assign it to norb = # of AOs
    is >> norb;
    cout << "Number of atomic orbitals: " << norb << endl;

    //Using norbs, create an norb x norb Matrix for the oei
    S.resize(norb, norb);
    T.resize(norb, norb);
    V.resize(norb, norb);
    core.resize(norb, norb);
    SOM.resize(norb, norb);

    //Create linear array for Two electron integral
    int M = (norb*(norb+1))/2;      //Number of elements in matrix ixj
    int N = (M*(M+1))/2;            //Number of elements in super matrix ijxkl
    TEI.resize(N);                  //So this will just be a linear array of size N to hold all elements of ijxkl

    //Create Fock Matrix of same dimension as core and oei matrices
    F_Guess.resize(norb,norb);

    is.close(); //Close input file
}


//Print out whatever matrix I assign to it (ixj matrix)
void HartreeFock::print_matrix(std::string mat_string, Matrix matrix)
{
    cout << endl;
    cout << mat_string;
    for(int i=0; i<matrix.rows(); i++) {
        for(int j=0; j<matrix.cols(); j++) {
            printf("%13.7f", matrix(i,j));
        }
        printf("\n");
    }
    cout << endl;
    return;
}
void HartreeFock::print_matrix(std::string mat_string, Vector vect)
{
    cout << endl;
    cout << mat_string;
    for(int i=0; i<vect.size(); i++) {
        printf("%13.7f\n", vect(i));
    }
    cout << endl;
    return;
}


//Read in and print out nuclear repulsion energy
void HartreeFock::read_enuc(const char *filename)
{
    std::ifstream nucl(filename);
    assert(nucl.good());
    nucl >> enuc;
    cout << endl;
    printf("Nuclear Repulsion Energy: %12.15f \n", enuc);
    nucl.close();
    return;
}


//Read in one electron integrals and store into respective matrices
Matrix HartreeFock::read_overlap(HartreeFock hf, const char *filename)
{
    //Open File here
    std::ifstream oei(filename);
    assert(oei.good());

    //Read in data
    int m;
    int n;
    while( oei >> m >> n >> hf.S(m-1,n-1) ) {
        hf.S(n-1,m-1) = hf.S(m-1,n-1);
    }
    oei.close(); //Close input file
    return hf.S;
}

Matrix HartreeFock::read_kinetic(HartreeFock hf, const char *filename)
{
    //Open File here
    std::ifstream oei(filename);
    assert(oei.good());

    //Read in data
    int m;
    int n;
    while( oei >> m >> n >> hf.T(m-1,n-1) ) {
        hf.T(n-1,m-1) = hf.T(m-1,n-1);
    }
   
    oei.close(); //Close input file
    return hf.T;
}

Matrix HartreeFock::read_potential(HartreeFock hf, const char *filename)
{
    //Open File here
    std::ifstream oei(filename);
    assert(oei.good());

    //Read in data
    int m;
    int n;
    while( oei >> m >> n >> hf.V(m-1,n-1) ) {
        hf.V(n-1,m-1) = hf.V(m-1,n-1);
    }
   
    oei.close(); //Close input file
    return hf.V;
}


//Build Core Hamiltonian from Kinetic Energy Integral and Nuclear Attraction Integral 
Matrix HartreeFock::build_core(HartreeFock hf)
{
    for(int i=0; i<core.rows(); i++) {
        for(int j=0; j<core.cols(); j++) {
            core(i,j) = hf.T(i,j) + hf.V(i,j);
        }
    }
    return core;
}


//Read in two-electron repulsion integral
Vector HartreeFock::read_tei(HartreeFock hf, const char *filename)
{
    //Open File here
    std::ifstream tei(filename);
    assert(tei.good());

    //Read in file
    int i, j, k, l, ij, kl, ijkl;
    double tei_val;        //Just need something to hold the value read in
    while( tei >> i >> j >> k >> l >> tei_val ) {

        i-=1;
        j-=1;
        k-=1;
        l-=1;
        
        if(i>j) ij = i*(i+1)/2 + j;
        else ij = j*(j+1)/2 + i;

        if(k>l) kl = k*(k+1)/2 + l;
        else kl = l*(l+1)/2 + k;
        
        if(ij>kl) ijkl = (ij*(ij+1)/2)+kl;
        else ijkl = (kl*(kl+1)/2)+ij;

        hf.TEI(ijkl) = tei_val;
    }
    tei.close(); //Close input file
    return hf.TEI;
}


//Build the orthogonalization matrix
Matrix HartreeFock::build_orthog(HartreeFock hf)
{
    //To diagonalize we need to solve for eigenvectors and eigenvalues
    Eigen::SelfAdjointEigenSolver<Matrix> solver(hf.S);
    Matrix evc = solver.eigenvectors();     //This is a matrix nxn
    Matrix evc_T = evc.transpose();         //This will stay nxn
    Matrix evl = solver.eigenvalues();      //This is a vector nx1

    //Take one over the squareroot of the eigenvalues
    for(int i=0; i<evl.size(); i++)
        evl(i) = 1/(sqrt(evl(i)));

    //Make sure the eigenvalues are a Diagonal Matrix
    Matrix evl_D = evl.asDiagonal();     //This should be nxn
    SOM = evc * evl_D * evc_T;

    return SOM;
}


//Build the Initial Guess Density
//
//First, form the Initial (guess) Fock Matrix
Matrix HartreeFock::build_fock_guess(HartreeFock hf)
{
    Matrix S_T = hf.S.transpose();      //Transpose of Symmetized Orthogonal Overlap Matrix (SOM)
    F_Guess = S_T * hf.core * hf.S;          // core_mat is Core Hamiltonian used as guess

    return F_Guess;
}

//Second, Diagonalize the Fock Matrix and transform its e-vectors into the og AO basis
Matrix HartreeFock::build_MO_coef(HartreeFock hf)
{
    //Diagonalize the Fock matrix
    Eigen::SelfAdjointEigenSolver<Matrix> solver(hf.F_Guess);
    Matrix C_p0 = solver.eigenvectors();   //The eigenvectors we will use in the transformation
    Matrix E_0 = solver.eigenvalues();     //The eigenvalules - E_0 matrix containing the initial orbital energies

    print_matrix("Eigenvectors (C' Matrix): \n", C_p0);
    print_matrix("Eigenvalues (Orbital energies): \n", E_0);

    //Transform the e-vectors into the original (non-orthogonal) AO basis
    MO_coef = hf.SOM * C_p0;
    
    return MO_coef;
}

//Third, build the Density Matrix using the occupied MOs
void HartreeFock::build_density(Matrix MO_coef, int elec_num) 
{   
    //We will have calculated total # of electrons in molecule.cc
    //Use that value to find the # of double occupied orbitals
    int m = elec_num/2;    // m is the number of doubly-occupied orbitals
    cout << "Number of doubley occupied orbitals: " << m << endl;

    return;
}


//Delete allocated and used memory
HartreeFock::~HartreeFock()
{
}

