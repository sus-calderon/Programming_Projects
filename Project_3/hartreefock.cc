// <hartreefock.cc>
// Written on July 16, 2019 by Susana Calderon
//

#include "hartreefock.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cassert>
#include <cmath>

//Include Eigen package for easy diagonalization
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;


//This will be to read in my integral. The enuc data will be it's own thing in the main program. (It's just one line and value.)
HartreeFock::HartreeFock(const char *filename)
{
    //Open File here
    std::ifstream is(filename);
    assert(is.good());
   
    //The next part is to go that last line and read in the first value which should tell me the total # of atoms
    //This should work if every line is the same length (checked with bash = 32 character total for each line)
    is.seekg(0,std::ios_base::end);      //So this should start at end of file...
    char ch = ' ';  //Init ch not equal to '\n'
    while(ch != '\n'){
        is.seekg(-2,std::ios_base::cur);   //2 steps back, meaning it doesn't check last character
        if((int)is.tellg() <= 0 ){  //If passed the start of the file, this is the start of the line
            is.seekg(0);
            break;
        }
        is.get(ch); //Checks the next character
    }

    //Now that we've accessed last line, we can grab first value and put it into natoms
    is >> norb;
    cout << "Number of atomic orbitals: " << norb << endl;

    //Now that I have norbs, I can create an norb X norb matrix
    //One electron integrals
    S = new double* [norb];
    for(int i=0; i<norb; i++)
        S[i] = new double[norb];

    T = new double* [norb];
    for(int i=0; i<norb; i++)
        T[i] = new double[norb];

    V = new double* [norb];
    for(int i=0; i<norb; i++)
        V[i] = new double[norb];

    core = new double* [norb];
    for(int i=0; i<norb; i++)
        core[i] = new double[norb];

    //Two electron integral
    //A linear array of size N
    int M = (norb*(norb+1))/2;      //Number of elements in matrix ixj
    int N = (M*(M+1))/2;            //Number of elements in super matrix ijxkl
    R = new double[N];            //So this will just be a linear array of size N to hold all elements of ijxkl

    is.close(); //Close input file
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

//Read in one electron integrals
void HartreeFock::read_oei(double** oei_mat, const char *filename)
{
    //Open File here
    std::ifstream oei(filename);
    assert(oei.good());

    //Read in data
    int m;
    int n;
    while( oei >> m >> n >> oei_mat[m-1][n-1] ) {
        oei_mat[n-1][m-1] = oei_mat[m-1][n-1];
    }
   
    oei.close(); //Close input file
    return;
}


//Print out one electron integrals we've read in
void HartreeFock::print_matrix(std::string mat_string, double** matrix)
{
    cout << endl;
    cout << mat_string;
    for(int i=0; i<norb; i++) {
        for(int j=0; j<norb; j++) {
            printf("%13.7f", matrix[i][j]);
        }
        printf("\n");
    }
    cout << endl;
    return;
}


// Build Core Function
void HartreeFock::build_core(double** t_mat, double** v_mat)
{
    for(int i=0; i<norb; i++) {
        for(int j=0; j<norb; j++) {
            core[i][j] = t_mat[i][j] + v_mat[i][j];
        }
    }

    return;
}


//Read in two-electron repulsion integral
void HartreeFock::read_tei(double* tei_ary, const char *filename)
{
    //Open File here
    std::ifstream tei(filename);
    assert(tei.good());

    //Read in file
    int i, j, k, l, ij, kl, ijkl;
    double tei_val;        //Just need something to hold the value read in
    while( tei>> i >> j >> k >> l >> tei_val ) {

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

        R[ijkl] = tei_val;
    }
    tei.close(); //Close input file

    return;
}


//Build the orthogonalization matrix
void HartreeFock::build_orthog(double** s_mat)
{
    //We Diagonalize the overlap matrix S
    //I can use the Eigen package to simplify what I write down in the code
    Matrix overlap(norb,norb);
    for(int i=0; i<norb; i++) {
        for(int j=0; j<norb; j++) {
            overlap(i,j) = s_mat[i][j];
        }
    }

    //To diagonalize we need to solve for eigenvectors and eigenvalues
    Eigen::SelfAdjointEigenSolver<Matrix> solver(overlap);
    Matrix EVC = solver.eigenvectors();   //This is a matrix nxn
    Matrix EVC_T = EVC.transpose();      //This will stay nxn
    Matrix EVL = solver.eigenvalues();    //This is a vector nx1

    //cout << endl;
    //cout << "S = Ls * D * Ls^(T) = " << endl << EVC * EVL.asDiagonal() * EVC_T << endl;

    //Take the squareroot of the eigenvalues
    for(int i=0; i<EVL.size(); i++)
        EVL(i) = sqrt(EVL(i));

    //Make sure the eigenvalues are a Diagonal Matrix
    Matrix EVL_D = EVL.asDiagonal();     //This should be nxn

    //Ask Kirk about Precision Issue
    cout << endl;
    cout << "S^1/2 = Ls * D^1/2 * Ls^(T) = " << endl << EVC * EVL_D * EVC_T << endl;

    return;
}


//Delete allocated and used memory
HartreeFock::~HartreeFock()
{
    for(int i=0; i<norb; i++)
        delete[] S[i];
    delete[] S;

    for(int i=0; i<norb; i++)
        delete[] T[i];
    delete[] T;

    for(int i=0; i<norb; i++)
        delete[] V[i];
    delete[] V;

    for(int i=0; i<norb; i++)
        delete[] core[i];
    delete[] core;

    delete[] R;
}


