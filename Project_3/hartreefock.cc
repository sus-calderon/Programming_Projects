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

    //Two electron integral
    //A linear array of size N
    int M = (norb*(norb+1))/2;      //Number of elements in matrix ixj
    int N = (M*(M+1))/2;            //Number of elements in super matrix ijxkl
    TEI = new double[N];            //So this will just be a linear array of size N to hold all elements of ijxkl

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


//Read in one electron integrals and store into respecitve matrices
Matrix HartreeFock::read_overlap(Matrix oei_mat, const char *filename)
{
    //Open File here
    std::ifstream oei(filename);
    assert(oei.good());

    //Read in data
    int m;
    int n;
    while( oei >> m >> n >> oei_mat(m-1,n-1) ) {
        oei_mat(n-1,m-1) = oei_mat(m-1,n-1);
    }
   
    oei.close(); //Close input file
    return oei_mat;
}

Matrix HartreeFock::read_kinetic(Matrix oei_mat, const char *filename)
{
    //Open File here
    std::ifstream oei(filename);
    assert(oei.good());

    //Read in data
    int m;
    int n;
    while( oei >> m >> n >> oei_mat(m-1,n-1) ) {
        oei_mat(n-1,m-1) = oei_mat(m-1,n-1);
    }
   
    oei.close(); //Close input file
    return oei_mat;
}

Matrix HartreeFock::read_potential(Matrix oei_mat, const char *filename)
{
    //Open File here
    std::ifstream oei(filename);
    assert(oei.good());

    //Read in data
    int m;
    int n;
    while( oei >> m >> n >> oei_mat(m-1,n-1) ) {
        oei_mat(n-1,m-1) = oei_mat(m-1,n-1);
    }
   
    oei.close(); //Close input file
    return oei_mat;
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


//Build Core Hamiltonian from Kinetic Energy Integral and Nuclear Attraction Integral 
Matrix HartreeFock::build_core(Matrix t_mat, Matrix v_mat)
{
    for(int i=0; i<core.rows(); i++) {
        for(int j=0; j<core.cols(); j++) {
            core(i,j) = t_mat(i,j) + v_mat(i,j);
        }
    }

    return core;
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

        TEI[ijkl] = tei_val;
    }
    tei.close(); //Close input file
    return;
}


//Build the orthogonalization matrix
void HartreeFock::build_orthog(Matrix s_mat)
{
    //To diagonalize we need to solve for eigenvectors and eigenvalues
    Eigen::SelfAdjointEigenSolver<Matrix> solver(s_mat);
    Matrix EVC = solver.eigenvectors();     //This is a matrix nxn
    Matrix EVC_T = EVC.transpose();         //This will stay nxn
    Matrix EVL = solver.eigenvalues();      //This is a vector nx1

    //Take one over the squareroot of the eigenvalues
    for(int i=0; i<EVL.size(); i++)
        EVL(i) = 1/(sqrt(EVL(i)));

    //Make sure the eigenvalues are a Diagonal Matrix
    Matrix EVL_D = EVL.asDiagonal();     //This should be nxn

    //cout << endl;
    //cout << "S^1/2 = Ls * D^1/2 * Ls^(T) = " << endl << EVC * EVL_D * EVC_T << endl;
    return;
}


//Build the Initial Guess Density
void HartreeFock::build_fock_guess(double** s_ortho, Matrix core_mat)
{

    return;
}


//Delete allocated and used memory
HartreeFock::~HartreeFock()
{
    delete[] TEI;
}

