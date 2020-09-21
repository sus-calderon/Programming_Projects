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

#define BIGNUM 1000     //So defining BIGNUM here allows it be of GLOBAL use
#define INDEX(i,j) (i>j) ? (ioff(i)+j) : (ioff(j)+i)
//#define INDEX(i,j) ((i>j) ? (((i)*((i)+1)/2)+ (j)) : (((j)*((j)+1/2)+(i)))

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
    //cout << "Number of atomic orbitals: " << norb << endl;

    //Using norbs, create an norb x norb Matrix for the oei
    S.resize(norb, norb);
    T.resize(norb, norb);
    V.resize(norb, norb);
    core.resize(norb, norb);
    SOM.resize(norb, norb);
    D.resize(norb,norb);

    //Create linear array for Two electron integral
    int M = (norb*(norb+1))/2;      //Number of elements in matrix ixj
    int N = (M*(M+1))/2;            //Number of elements in super matrix ijxkl
    TEI.resize(N);                  //So this will just be a linear array of size N to hold all elements of ijxkl
    TEI_MO.resize(N);

    ioff.resize(BIGNUM);

    //Create Fock Matrix of same dimension as core and oei matrices
    F.resize(norb,norb);

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

void HartreeFock::print_vector(std::string mat_string, Vector vect)
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
double HartreeFock::read_enuc(const char *filename)
{
    std::ifstream nucl(filename);
    assert(nucl.good());
    nucl >> enuc;
    cout << endl;
    //printf("Nuclear Repulsion Energy: %12.15f \n", enuc);
    nucl.close();
    return enuc;
}


//Read in one electron integrals and store into respective matrices
int HartreeFock::read_overlap(HartreeFock& hf, const char *filename)
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

    //print_matrix("Overlap Integral Matrix (s): \n", hf.S);

    return 0;
}

int HartreeFock::read_kinetic(HartreeFock& hf, const char *filename)
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

    //hf.print_matrix("Kinetic Energy Integral Matrix (t): \n", hf.T);

    return 0;
}

int HartreeFock::read_potential(HartreeFock& hf, const char *filename)
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

    //hf.print_matrix("Nuclear Attraction Integral Matrix (v): \n", hf.V);

    return 0;
}


//Build Core Hamiltonian from Kinetic Energy Integral and Nuclear Attraction Integral 
int HartreeFock::build_core(HartreeFock& hf)
{
    for(int i=0; i<hf.core.rows(); i++) {
        for(int j=0; j<hf.core.cols(); j++) {
            hf.core(i,j) = hf.T(i,j) + hf.V(i,j);
        }
    }

    //hf.print_matrix("Core Hamiltonian Matrix (h): \n", hf.core);

    return 0;
}


//Read in two-electron repulsion integral
int HartreeFock::read_tei(HartreeFock& hf, const char *filename)
{
    //Open File here
    std::ifstream tei(filename);
    assert(tei.good());

    //Read in file
    ioff(0) = 0;
    for(int n=1; n<1000; n++) {
       ioff(n) = ioff(n-1) + n;
    }

    int i, j, k, l, ij, kl, ijkl;
    double tei_val;        //Just need something to hold the value read in
    while( tei >> i >> j >> k >> l >> tei_val ) {
        i-=1;
        j-=1;
        k-=1;
        l-=1;
        
        //ij = (i>j) ? (i*(i+1)/2 + j) : (j*(j+1)/2 + i);
        //This is a shortened version of the if else statement below
        // using conditional (ternary) operators
        //(condition) ? (if_true) : (if_false);
        //if(i>j) ij = i*(i+1)/2 + j;
        //else ij = j*(j+1)/2 + i;

        ij = (i>j) ? (ioff(i) + j) : (ioff(j) + i);
        kl = (k>l) ? (ioff(k) + l) : (ioff(l) + k);
        ijkl = (ij>kl) ? (ioff(ij) + kl) : (ioff(kl) + ij); 

        hf.TEI(ijkl) = tei_val;

    }
    tei.close(); //Close input file

    //hf.print_vector("TEI array: \n", hf.TEI);

    return 0;
}


//Build the orthogonalization matrix
int HartreeFock::build_orthog(HartreeFock& hf)
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
    hf.SOM = evc * evl_D * evc_T;

    //hf.print_matrix("Symmetric Orthogonalization Matrix (S^1/2): \n", hf.SOM); 

    return 0;
}


//This is Step 7, to which we return to if the difference in consecutive SCF energy and the 
//root-mean-squared difference in consecutive densities do not fall below the prescribed thresholds.
//
// Step 7: Compute the new Fock matrix (F) for the SCF procedure
int HartreeFock::update_Fock(HartreeFock& hf)
{
    int i, j, k, l, ij, kl, ijkl, ik, jl, ikjl;
    hf.F = hf.core;                    //So the new Fock matrix comes from adding core_H to Density*TEI
    for(i=0; i<hf.F.rows(); i++) {
        for(j=0; j<hf.F.rows(); j++) {
            for(k=0; k<hf.F.rows(); k++) {
                for(l=0; l<hf.F.rows(); l++) {
                    ij = INDEX(i,j);
                    //if(i>j) ij = i*(i+1)/2 + j;
                    //else ij = j*(j+1)/2 + i;

                    kl = INDEX(k,l);
                    ijkl = INDEX(ij,kl);
                    ik = INDEX(i,k);
                    jl = INDEX(j,l);
                    ikjl = INDEX(ik,jl);

                    hf.F(i,j) += hf.D(k,l) * (2.0 * hf.TEI(ijkl) - hf.TEI(ikjl));
                }
            }
        } 
    }

    if(hf.iter == 1){
        //hf.print_matrix("Fock Matrix (F): \n", hf.F);
    }
    
    return 0;
}
//
// Step 8: Build the New Density Matrix
int HartreeFock::build_density(HartreeFock& hf, int elec_num) 
{   
    //Transform Fock matrix (F -> F')
    hf.F_p = hf.SOM.transpose() * hf.F * hf.SOM;            //Now new Fock matrix is used as guess 

    //Print F' on first iteration
    if(hf.iter == 0){
        //hf.print_matrix("Initial Fock Matrix (F'): \n", hf.F_p); 
    }

    //Diagonalize the F'
    Eigen::SelfAdjointEigenSolver<Matrix> solver(hf.F_p);
    Matrix C_p = solver.eigenvectors();   //The eigenvectors we will use in the transformation
    //Matrix E = solver.eigenvalues();     //The eigenvalules - E_0 matrix containing the initial orbital energies
    hf.E_p = solver.eigenvalues();

    //print_matrix("Eigenvectors (C' Matrix): \n", C_p0);
    //print_vector("Eigenvalues (Orbital energies): \n", E_0);

    //Transform the e-vectors into the original (non-orthogonal) AO basis
    hf.C = hf.SOM * C_p;
    
    //Print C Matrix on first iteration
    if(hf.iter == 0){
        //hf.print_matrix("Initial Coefficient Matrix (C): \n", hf.C);
    }

    //Build Density
    //We will have calculated total # of electrons in molecule.cc
    int occ = elec_num/2;             // occ is the number of doubly-occupied orbitals
    Matrix C_d = hf.C.block(0,0,hf.C.rows(),occ);
    hf.D = C_d * C_d.transpose();      //C_do is 7x5 and C_do_T is 5x7 so my resulting matrix is 7x7

    //Print density on first iteration
    if(hf.iter == 0){
        //hf.print_matrix("Initial Density Matrix (D): \n", hf.D);
    }

    return 0;
}
//
// Step 9: Compute the New SCF Energy
int HartreeFock::compute_SCF(HartreeFock& hf)
{
    hf.SCF = 0.0;
    for(int i=0; i<hf.D.rows(); i++) {
        for(int j=0; j<hf.D.cols(); j++) {
            hf.SCF += hf.D(i,j) * (hf.core(i,j) + hf.F(i,j));
        }
    }
    hf.tot_E = hf.SCF + hf.enuc;
    return 0;
}

// MP2 CODE
// Step 3: Transform the 2e integrals from AO to MO basis
int HartreeFock::transform_AO_2_MO(HartreeFock& hf)
{
    int ijkl;
    int pq, rs, pqrs;

    //Noddy code:
    for(int i=0; i<norb; i++)
    {
        for(int j=0; j<=i; j++)
        {
            for(int k=0; k<=i; k++)
            {
                for(int l=0; l<=(i==k ? j:k); l++, ijkl++)
                {
                    for(int p=0; p<norb; p++)
                    {
                        for(int q=0; q<norb; q++)
                        {
                            pq = INDEX(p,q);
                            for(int r=0; r<norb; r++)
                            {
                                for(int s=0; s<norb; s++) {
                                    rs=INDEX(r,s);
                                    pqrs=INDEX(pq,rs);
                                    //no use of hf.TEI_MO since it's something I barely created in this function
                                    TEI_MO(ijkl) += hf.C(p,i) * hf.C(q,j) * hf.C(r,k) * hf.C(s,l) * hf.TEI(pqrs);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return 0;
}

int HartreeFock::MP2_calc(HartreeFock& hf, int elec_num)
{
    //do I have to wonder about it being float or double?
    double Emp2 = 0.0;
    int ndocc = elec_num/2;
    //int nao = norb;  
    int ia, ja, jb, ib, iajb, ibja;
    //OK think about # of doubly occupied vs unocc/virtual orbitals(these come after)
    //eps are the orbital energies I solved for as eigenvalues
    for(int i=0; i < ndocc; i++) {
        for(int a=ndocc; a<norb; a++) {
            ia = INDEX(i,a);
            for(int j=0; j<ndocc; j++) {
                ja = INDEX(j,a);
                for(int b=ndocc; b<norb; b++) {
                    jb = INDEX(j,b);
                    ib = INDEX(i,b);
                    iajb = INDEX(ia, jb);
                    ibja = INDEX(ib, ja);
                    Emp2 += (hf.TEI_MO[iajb] * (2*hf.TEI_MO[iajb] - hf.TEI_MO[ibja]))/(hf.E_p(i) + hf.E_p(j) - hf.E_p(a) - hf.E_p(b));
                    //eps is epsilon  which are orbital energies that were solved for using eigen solver
                }
            }
        }
    }
    cout << "MP2 Energy: " << Emp2 << std::endl;
    return 0;
}

//CCSSD Functions



//Delete allocated and used memory
HartreeFock::~HartreeFock()
{
}

