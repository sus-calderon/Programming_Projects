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
HartreeFock::HartreeFock(const char *filename, int e_c)
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
    nmo = 2*norb;   //# of spin MOs
    no = e_c;     //# of occupied spin orbitals
    nv = nmo-no;    //# of virtual orbitals (was previously e_uo)

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

    tmp.resize((norb*(norb+1)/2),(norb*(norb+1)/2));
    X.resize(norb,norb);
    Y.resize(norb,norb);

    ioff.resize(BIGNUM);

    //Create Fock Matrix of same dimension as core and oei matrices
    F.resize(norb,norb);

    //Begin for-loop for creating 4D array (and assign its dimensions)
    spinTEI = new double***[nmo];
    for(int i=0; i<nmo; i++) {
        spinTEI[i] = new double**[nmo];
        for(int j=0; j<nmo; j++) {
            spinTEI[i][j] = new double*[nmo];
            for(int k=0; k<nmo; k++) {
                spinTEI[i][j][k] = new double[nmo];
                for(int l=0; l<nmo; l++) {
                    spinTEI[i][j][k][l] = 0.0; //Initialize matrix
                }
            }
        }
    }

    int e_uo = nmo - e_c; //unoccupied electrons that would be virtual
    T2 = new double***[nmo];
    for(int i=0; i<nmo; i++) {
        T2[i] = new double**[nmo];
        for(int j=0; j<nmo; j++) {
            T2[i][j] = new double*[nmo];
            for(int a=0; a<nmo; a++) {
                T2[i][j][a] = new double[nmo];
                for(int b=0; b<nmo; b++) {
                    T2[i][j][a][b] = 0.0;
                }
            }
        }
    }

    //create effective 2-particle excitation operator tau and tau_p
    tau = new double***[nmo];
    for(int i=0; i<nmo; i++) {
        tau[i] = new double**[nmo];
        for(int j=0; j<nmo; j++) {
            tau[i][j] = new double *[nmo];
            for(int a=0; a<nmo; a++) {
                tau[i][j][a] = new double[nmo];
                for(int b=0; b<nmo; b++) {
                    tau[i][j][a][b] = 0.0;
                }
            }
        }
    }

    tau_t = new double***[nmo];
    for(int i=0; i<nmo; i++) {
        tau_t[i] = new double**[nmo];
        for(int j=0; j<nmo; j++) {
            tau_t[i][j] = new double*[nmo];
            for(int a=0; a<nmo; a++) {
                tau_t[i][j][a] = new double[nmo];
                for(int b=0; b<nmo; b++) {
                    tau_t[i][j][a][b] = 0.0;
                }
            }
        }
    }

    T1.resize(nmo, nmo);
    T1.setZero();
    
    D_ia.resize(nmo, nmo);
    D_ia.setZero();
    
    //F Intermeidates
    F_ae.resize(nmo, nmo);
    F_mi.resize(nmo, nmo);
    F_me.resize(nmo, nmo);

    F_ae.setZero();
    F_mi.setZero();
    F_me.setZero();

    //W Intermediates
    W_mnij = new double***[nmo];
    for(int i=0; i<nmo; i++) {
        W_mnij[i] = new double**[nmo];
        for(int j=0; j<nmo; j++) {
            W_mnij[i][j] = new double*[nmo];
            for(int a=0; a<nmo; a++) {
                W_mnij[i][j][a] = new double[nmo];
                for(int b=0; b<nmo; b++) {
                    W_mnij[i][j][a][b] = 0.0;
                }
            }
        }
    }

    W_abef = new double***[nmo];
    for(int i=0; i<nmo; i++) {
        W_abef[i] = new double**[nmo];
        for(int j=0; j<nmo; j++) {
            W_abef[i][j] = new double*[nmo];
            for(int a=0; a<nmo; a++) {
                W_abef[i][j][a] = new double[nmo];
                for(int b=0; b<nmo; b++) {
                    W_abef[i][j][a][b] = 0.0;
                }
            }
        }
    }

    W_mbej = new double***[nmo];
    for(int i=0; i<nmo; i++) {
        W_mbej[i] = new double**[nmo];
        for(int j=0; j<nmo; j++) {
            W_mbej[i][j] = new double*[nmo];
            for(int a=0; a<nmo; a++) {
                W_mbej[i][j][a] = new double[nmo];
                for(int b=0; b<nmo; b++) {
                    W_mbej[i][j][a][b] = 0.0;
                }
            }
        }
    }

    T2 = new double***[nmo];
    for(int i=0; i<nmo; i++) {
        T2[i] = new double**[nmo];
        for(int j=0; j<nmo; j++) {
            T2[i][j] = new double*[nmo];
            for(int a=0; a<nmo; a++) {
                T2[i][j][a] = new double[nmo];
                for(int b=0; b<nmo; b++) {
                    T2[i][j][a][b] = 0.0;
                }
            }
        }
    }

    D_ijab = new double***[nmo];
    for(int i=0; i<nmo; i++) {
        D_ijab[i] = new double**[nmo];
        for(int j=0; j<nmo; j++) {
            D_ijab[i][j] = new double*[nmo];
            for(int a=0; a<nmo; a++) {
                D_ijab[i][j][a] = new double[nmo];
                for(int b=0; b<nmo; b++) {
                    D_ijab[i][j][a][b] = 0.0;
                }
            }
        }
    }

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
    //cout << norb << endl;
    //cout << e_count << endl;
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
    
    Eigen::SelfAdjointEigenSolver<Matrix> solver(hf.F);
    hf.F_eval = solver.eigenvalues();

    return 0;
}

// Step 8: Build the New Density Matrix
int HartreeFock::build_density(HartreeFock& hf, int elec_num) 
{   
    //Transform Fock matrix (F -> F')
    hf.F_p = hf.SOM.transpose() * hf.F * hf.SOM;            //Now new Fock matrix is used as guess 

    //Print F' on first iteration
    //if(hf.iter == 0){
        //hf.print_matrix("Initial Fock Matrix (F'): \n", hf.F_p); 
    //}

    //Diagonalize the F'
    Eigen::SelfAdjointEigenSolver<Matrix> solver(hf.F_p);
    Matrix C_p = solver.eigenvectors();   //The eigenvectors we will use in the transformation
    hf.E_p = solver.eigenvalues();

    //Transform the e-vectors into the original (non-orthogonal) AO basis
    hf.C = hf.SOM * C_p;
    
    //Print C Matrix on first iteration
    //if(hf.iter == 0){
        //hf.print_matrix("Initial Coefficient Matrix (C): \n", hf.C);
    //}

    //Build Density
    int occ = elec_num/2;             // occ is the number of doubly-occupied orbitals
    Matrix C_d = hf.C.block(0,0,hf.C.rows(),occ);
    hf.D = C_d * C_d.transpose();      //C_do is 7x5 and C_do_T is 5x7 so my resulting matrix is 7x7

    //Print density on first iteration
    //if(hf.iter == 0){
        //hf.print_matrix("Initial Density Matrix (D): \n", hf.D);
    //}

    return 0;
}

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

//Smart Algorithm
//This transofrmation from AO to MO basis reduces the N^8 computational cost down to N^5
int HartreeFock::smart_transform_AO_2_MO(HartreeFock& hf)
{
    //hf.tmp.setZero((norb*(norb+1)/2), (norb*(norb+1)/2)); //this matrix holds in values from X and Y, where X and Y are each TEI and when combined needs to account for their dimensions

    int i,j,k,l, ij,kl,ijkl,klij;

    for(i=0, ij=0; i<norb; i++) 
        for(j=0; j<=i; j++,ij++) {
            for(k=0, kl=0; k<norb; k++) 
                for(l=0; l<=k; l++,kl++) {
                    ijkl=INDEX(ij,kl);
                    hf.X(k,l)=hf.X(l,k)=hf.TEI(ijkl);
                }
            //Matrix Y is already set to zero, with X containing TEI_AO values
            //Now do some Matrix multiplication
            hf.Y.setZero(norb,norb);
            hf.Y = hf.C.transpose()*hf.X;
            hf.X.setZero(norb,norb);
            hf.X = hf.Y*hf.C;

            for(k=0,kl=0; k<norb; k++)
                for(l=0; l<=k; l++, kl++)
                    hf.tmp(kl,ij) = hf.X(k,l);
        }

    hf.TEI.setZero((norb*(norb+1)/2)*((norb*(norb+1)/2)+1)/2);
    for(k=0, kl=0; k<norb; k++)
        for(l=0; l<=k; l++, kl++) {
            hf.X.setZero(norb,norb);
            hf.Y.setZero(norb,norb);
            for(i=0, ij=0; i<norb; i++)
                for(j=0; j<=i; j++,ij++)
                    hf.X(i,j) = hf.X(j,i) = hf.tmp(kl,ij);
            hf.Y.setZero(norb,norb);
            hf.Y = hf.C.transpose()*hf.X;
            hf.X.setZero(norb,norb);
            hf.X = hf.Y*hf.C;
            for(i=0, ij=0; i<norb; i++)
                for(j=0; j<=i; j++, ij++) {
                    klij = INDEX(kl,ij);
                    hf.TEI(klij) = hf.X(i,j);
                }
        }

    return 0;
}


int HartreeFock::smart_MP2_calc(HartreeFock& hf, int elec_num)
{
    //do I have to wonder about it being float or double?
    hf.Emp2 = 0.0;
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
                    hf.Emp2 += (hf.TEI[iajb] * (2*hf.TEI[iajb] - hf.TEI[ibja]))/(hf.E_p(i) + hf.E_p(j) - hf.E_p(a) - hf.E_p(b));
                    //eps is epsilon  which are orbital energies that were solved for using eigen solver
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
                    Emp2 += (hf.spinTEI[i][a][j][b] * (2*hf.spinTEI[i][a][j][b] - hf.spinTEI[i][b][j][a]))/(hf.E_p(i) + hf.E_p(j) - hf.E_p(a) - hf.E_p(b));
                    //eps is epsilon  which are orbital energies that were solved for using eigen solver
                }
            }
        }
    }
    //cout << "MP2 Energy: " << Emp2 << std::endl;
    return 0;
}

//CCSD Functions
void HartreeFock::transform_to_spin(HartreeFock& hf)
{
    int pr, ps, qr, qs, prqs, psqr;
    double val1, val2;
    //again, in order to access all orbitals, whether they are spin up or spin down, I multiply norb by 2
    for(int p=0; p<2*norb; p++) {
        for(int q=0; q<2*norb; q++) {
            for(int r=0; r<2*norb; r++) {
                for(int s=0; s<2*norb; s++) {
                    //Start Indexing
                    pr = INDEX(p/2, r/2);
                    qs = INDEX(q/2, s/2);
                    prqs = INDEX(pr, qs);
                    //val1 = value of integral <pq|rs> when p=r and q=s
                    val1 = (hf.TEI[prqs])*(p%2==r%2)*(q%2==s%2);
                    //Next Indexing for psqr
                    ps = INDEX(p/2, s/2);
                    qr = INDEX(q/2, r/2);
                    psqr = INDEX(ps, qr);
                    //val2 = <pq|sr> when p=s and q=r
                    val2 = (hf.TEI[psqr])*(p%2==s%2)*(q%2==r%2);
                    //Now combine to get <pq||rs> = antisym MO spin orbital basis
                    hf.spinTEI[p][q][r][s] = val1 - val2;
                }
            }
        }
    }
    return;
}

//Create spin-orbital Fock matrix
Matrix HartreeFock::create_spinFock(HartreeFock& hf, int elec_num) 
{
    hf.spinFock.setZero(2*norb,2*norb);
    Matrix core_MO(norb,norb);
    core_MO = hf.C.transpose()*hf.core*hf.C;      //Transform the core_H into the MO basis

    //Begin for loops
    for(int p=0; p<2*norb; p++) {
        for(int q=0; q<2*norb; q++) {
           hf.spinFock(p,q) = core_MO(p/2, q/2)*(p%2 == q%2);
           for(int m=0; m<elec_num; m++) {
               hf.spinFock(p,q) += hf.spinTEI[p][m][q][m];
           }
           if(std::abs(hf.spinFock(p,q)) < std::pow(10,-7)) {
               hf.spinFock(p,q) = 0.0;
           }
        }
    }
    cout<<"The Fock matrix in the spin orbital basis:"<<std::endl<<hf.spinFock<<std::endl;
    return hf.spinFock ;
}

//Build inital-guess cluster amplitudes
void HartreeFock::build_cluster_amp(HartreeFock& hf, int e_n)
{
    //These are your doubles
    // <ij||ab> = <ij|ab> - <ij|ba>
    int e_uo = (2*norb) - e_n; //unoccupied electrons that would be virtual
    for(int i=0; i<e_n; i++) {
        for(int j=0; j<e_n; j++) {
            for(int a=0; a<e_uo; a++) {
                for(int b=0; b<e_uo; b++) {
                    hf.T2[i][j][a][b] = (hf.spinTEI[i][j][a+e_n][b+e_n])/(hf.spinFock(i,i) + hf.spinFock(j,j) - hf.spinFock(a+e_n,a+e_n) - hf.spinFock(b+e_n,b+e_n));
                    //I have to add e_n to a and b in order to access electrons within the virtual orbitals (otherwise I'm calling for occupied)
                }
            }
        }
    }
    //double check that the MP2 energy is the same as in project 4
    hf.spinEmp2 = 0.0;
    for(int i=0; i<e_n; i++) {
        for(int j=0; j<e_n; j++) {
            for(int a=0; a<e_uo; a++) {
                for(int b=0; b<e_uo; b++) {
                    hf.spinEmp2 += (hf.T2[i][j][a][b])*(hf.spinTEI[i][j][a+e_n][b+e_n]);
                }
            }
        }
    }
    hf.spinEmp2*=0.25;
    cout << "The MP2 energy correction using the spin orbital basis: " << hf.spinEmp2 << "\n" << endl;

    return;
}

void HartreeFock::build_tau(HartreeFock& hf)
{
    //Tau has ij and ab indices so another 4D array must be built
    int e_c = e_count;
    int e_uo = (2*norb) - e_c;
    int i,j,a,b;
    for (i=0; i<e_c; i++) {
        for(j=0; j<e_c; j++) {
            for(a=0; a<e_uo; a++) {
                for(b=0; b<e_uo; b++) {
                    hf.tau[i][j][a][b] = hf.T2[i][j][a][b] + hf.T1(i,a)*hf.T1(j,b) - hf.T1(i,b)*hf.T1(j,a);
                    hf.tau_t[i][j][a][b] = hf.T2[i][j][a][b] + 0.5*(hf.T1(i, a)*hf.T1(j,b) - hf.T1(i,b)*hf.T1(j,a));
                }
            }
        }
    }
    return;
}

//Calculate Coupled Cluster Intermediates
void HartreeFock::ccF_intermediates(HartreeFock& hf, int e_c)
{
    //int nmo = 2*norb;   //# of spin MOs
    //int no = 2*e_c;     //# of occupied orbitals
    //int nv = nmo-no;    //# of virtual orbitals (was previously e_uo)

    double sum1 = 0.0;
    double sum2 = 0.0;
    double sum3 = 0.0;

    for(int a=0; a<nv; a++) {
        for(int e=0; e<nv; e++) {
            for(int m=0; m<no; m++) {
                sum1 += hf.spinFock(m,no+e)*hf.T1(m,a);
                
                for(int f=0; f<nv; f++) {
                    sum2 += hf.T1(m,f)*hf.spinTEI[m][no+a][no+f][no+e];
                    
                    for(int n=0; n<no; n++) {
                        sum3 += hf.tau_t[m][n][a][f] * hf.spinTEI[m][n][no+e][no+f];
                    }
                }
            }
            hf.F_ae(a,e) = (1-(a==e))*hf.spinFock(no+a,no+e) - 0.5*sum1 + sum2 - 0.5*sum3;
        }
    }
    
    sum1 = 0.0;
    sum2 = 0.0;
    sum3 = 0.0;

    for(int m=0; m<no; m++) {
        for(int i=0; i<no; i++) {
            for(int e=0; e<nv; e++) {
                sum1 += hf.spinFock(m,no+e)*T1(i,e);

                for(int n=0; n<no; n++) {
                    sum2 += hf.T1(n, e)*hf.spinTEI[m][n][i][no+e];

                    for(int f=0; f<nv; f++) {
                        sum3 += hf.tau_t[i][n][e][f] * hf.spinTEI[m][n][no+e][no+f];
                    }
                }
            }
            hf.F_mi(m,i) = (1-(m==i))*hf.spinFock(m,i) + 0.5*sum1 + sum2 + 0.5*sum3;
        }
    }


    sum1 = 0.0;

    for(int m=0; m<no; m++) {
        for(int e=0; e<nv; e++) {
            for(int n=0; n<no; n++) {
                for(int f=0; f<nv; f++) {
                    sum1 += hf.T1(n,f)*hf.spinTEI[m][n][no+e][no+f];
                }
            }
            hf.F_me(m,e) = hf.spinFock(m,no+e) + sum1;
        }
    }

    return;
}

void HartreeFock::ccW_intermediates(HartreeFock& hf, int e_c)
{
    //int nmo = 2*norb;   //# of spin MOs
    //int no = 2*e_c;     //# of occupied electrons
    //int nv = nmo-no;    //# of virtual electrons (was previously e_uo)

    //Wmnij
    double sum1 = 0.0;
    double sum2 = 0.0;
    for(int m=0; m<no; m++) {
        for(int n=0; n<no; n++) {
            for(int i=0; i<no; i++) {
                for(int j=0; j<no; j++) {
                    for(int e=0; e<nv; e++) {
                        //So before this sum, there is a Permutation operator for ij
                        //I'm guessing this switches ij to ji
                        sum1 += hf.T1(j,e)*hf.spinTEI[m][n][i][e+no];
                        sum1 -= hf.T1(i,e)*hf.spinTEI[m][n][j][e+no];
                        for(int f=0; f<nv; f++) {
                            sum2 += hf.tau[i][j][e][f]*hf.spinTEI[m][n][no+e][no+f];
                        }
                    }
                    hf.W_mnij[m][n][i][j] = hf.spinTEI[m][n][i][j] + sum1 + 0.25*sum2;
                }
            }
        }
    }

    //Wabef
    sum1 = 0.0;
    sum2 = 0.0;
    for(int a=0; a<nv; a++) {
        for(int b=0; b<nv; b++) {
            for(int e=0; e<nv; e++) {
                for(int f=0; f<nv; f++) {
                    for(int m=0; m<no; m++) {
                        sum1 += hf.T1(m,b)*hf.spinTEI[a+no][m][e+no][f+no];
                        sum1 -= hf.T1(m,a)*hf.spinTEI[b+no][m][e+no][f+no];
                        for(int n=0; n<no; n++) {
                            sum2 += hf.tau[m][n][a][b]*hf.spinTEI[m][n][no+e][no+f];
                        }
                    }
                    hf.W_abef[a][b][e][f] = hf.spinTEI[a+no][b+no][e+no][f+no] - sum1 + 0.25*sum2;
                }
            }
        }
    }

    //Wmbej
    sum1 = 0.0;
    sum2 = 0.0;
    double sum3 = 0.0;
    for(int m=0; m<no; m++) {
        for(int b=0; b<nv; b++) {
            for(int e=0; e<nv; e++) {
                for(int j=0; j<no; j++) {
                    for(int f=0; f<nv; f++) {
                        sum1 += hf.T1(j,f)*hf.spinTEI[m][b+no][e+no][f+no];
                        for(int n=0; n<no; n++) {
                            sum2 += hf.T1(n,b)*hf.spinTEI[m][n][no+e][j];
                            sum3 += (0.5*hf.T2[j][n][f][b] + hf.T1(j,f)*hf.T1(n,b))*hf.spinTEI[m][n][e+no][f+no];
                        }
                    }
                    hf.W_mbej[m][b][e][j] = hf.spinTEI[m][b+no][e+no][j] + sum1 - sum2 - sum3;
                }
            }
        }
    }

    return;
}

//Updating CC Terms
void HartreeFock::update_t_ia(HartreeFock& hf)
{
    //denominator arrays: D_ia = f_ii - f_aa    (diagonal terms occupied - virtual)
    // T1 Equation: t_ia*D_ia
    for(int i=0; i<no; i++) {
        for(int a=0; a<nv; a++) {
        double sum1 = 0.0;
        double sum2 = 0.0;
        double sum3 = 0.0;
        double sum4 = 0.0;
        double sum5 = 0.0;
        double sum6 = 0.0;
            for(int e=0; e<nv; e++) {
                sum1 += T1(i,e)*F_ae(a,e);
            }
            for(int m=0; m<no; m++) {
                    sum2 += T1(m,a)*F_mi(m,i);
            }
             
            for(int e=0; e<nv; e++) {
                for(int m=0; m<no; m++) {
                    sum3 += T2[i][m][a][e]*F_me(m,e);

                    for(int n=0; n<no; n++) {
                        sum6 += T2[n][m][a][e]*hf.spinTEI[n][m][e+no][i];
                    }
                }
            }
            for(int f=0; f<nv; f++) {
                for(int n=0; n<no; n++) {
                    sum4 += T1(n,f)*hf.spinTEI[n][a+no][i][f+no];
                }
            }
            for(int m=0; m<no; m++) {
                for(int e=0; e<nv; e++) {
                    for(int f=0; f<nv; f++) {
                        sum5 += T2[i][m][e][f]*hf.spinTEI[m][a+no][e+no][f+no];
                    }
                }
            }
            D_ia(i,a) = hf.spinFock(i,i) - hf.spinFock(a,a);
            T1(i,a) = hf.spinFock(i,a) + sum1 - sum2 + sum3 - sum4 - 0.5*sum5 - 0.5*sum6;  
            T1(i,a) /= D_ia(i,a);
        }
    }
    return;
}

void HartreeFock::update_t_ijab(HartreeFock& hf)
{
    // D_ijab = f_ii + f_jj - f_aa - f_bb   (for use with W intermediates)
    // T2 Equation: t_ijab*D_ijab
    for(int i=0; i<no; i++) {
        for(int j=0; j<no; j++) {
            for(int a=0; a<nv; a++) {
                for(int b=0; b<nv; b++) {
                    double sum1 = 0.0;
                    double sum2 = 0.0;
                    double sum3 = 0.0;
                    double sum4 = 0.0;
                    double sum5 = 0.0;
                    double sum6 = 0.0;
                    double sum7 = 0.0;
                    double sum8 = 0.0;
                    double sum9 = 0.0;
                    double sum10 = 0.0;
                    
                    for(int e=0; e<nv; e++) {
                        //P_(ab)
                        sum1 += T2[i][j][a][e]*F_ae(b,e);
                        sum1 -= T2[i][j][b][e]*F_ae(a,e);

                        for(int m=0; m<no; m++) {
                            //P_(ab)
                            sum2 += T2[i][j][a][e]*T1(m,b)*F_ae(b,e);
                            sum2 -= T2[i][j][b][e]*T1(m,a)*F_ae(a,e);
                        }
                    }

                    for(int m=0; m<no; m++) { 
                        //P_(ij)
                        sum3 += T2[i][m][a][b]*F_me(m,j);
                        sum3 -= T2[j][m][a][b]*F_me(m,i);

                        for(int e=0; e<nv; e++) { 
                            //P_(ij)
                            sum4 += T2[i][m][a][b]*(T1(j,e)*F_me(m,e));
                            sum4 -= T2[j][m][a][b]*(T1(i,e)*F_me(m,e));
                        }
                    }
                    
                    for(int m=0; m<no; m++) {
                        for(int n=0; n<no; n++) {
                            for(int e=0; e<nv; e++) {
                                for(int f=0; f<nv; f++) {
                                    //tau & W
                                    sum5 += tau[m][n][a][b]*hf.W_mnij[m][n][i][j];
                                    sum6 += tau[i][j][e][f]*hf.W_abef[a][b][e][f];
                                }
                            }
                        }
                    }

                    for(int m=0; m<no; m++) {
                        for(int e=0; e<nv; e++) {
                            //P_(ij) then P_(ab) on both i and j parts
                            sum7 += T2[i][m][a][e]*hf.W_mbej[m][b][e][j];
                            sum7 -= T2[j][m][a][e]*hf.W_mbej[m][b][e][i];
                            
                            //Now do P_(ab)
                            sum7 -= T2[i][m][b][e]*hf.W_mbej[m][a][e][j];
                            sum7 += T2[j][m][b][e]*hf.W_mbej[m][a][e][i];
                        }
                    }

                    for(int m=0; m<no; m++) {
                        for(int e=0; e<nv; e++) {
                            //P_(ij)
                            sum8 += T1(i,e)*T1(m,a)*(hf.spinTEI[m][b+no][e+no][j]);
                            sum8 -= T1(j,e)*T1(m,a)*(hf.spinTEI[m][b+no][e+no][i]);
                            
                            //P_(ab)
                            sum8 -= T1(i,e)*T1(m,b)*(hf.spinTEI[m][a+no][e+no][j]);
                            sum8 += T1(j,e)*T1(m,b)*(hf.spinTEI[m][a+no][e+no][i]);
                        }
                    }

                    for(int e=0; e<nv; e++) {
                        for(int m=0; m<no; m++) {
                            //P_(ij)
                            sum9 += T1(i,e)*hf.spinTEI[a+no][b+no][e+no][j];
                            sum9 -= T1(j,e)*hf.spinTEI[a+no][b+no][e+no][i];
                        }
                    }

                    for(int m=0; m<no; m++) {
                        //P(ab)
                        sum10 += T1(m,a)*hf.spinTEI[m][b][i][j];
                        sum10 -= T1(m,b)*hf.spinTEI[m][a][i][j];
                    }
                    //Use unsupported eigen package for 4D array
                    D_ijab[i][j][a][b] = hf.spinFock(i,i) + hf.spinFock(j,j) - hf.spinFock(a,a) - hf.spinFock(b,b);
                    T2[i][j][a][b] = hf.spinTEI[i][j][a+no][b+no] + (sum1 - 0.5*sum2) - (sum3 + 0.5*sum4) + 0.5*sum5 + 0.5*sum6 + (sum7 - sum8) + sum9 - sum10;
                    T2[i][j][a][b] /= D_ijab[i][j][a][b];
                }
            }
        }
    }

    return;
}

//Calculate CC Energy to check for convergence
//void HartreeFock::cc_E_convergence()
//{
//    double e_cc = 0.0;
//    double sum1 = 0.0;
//    double sum2 = 0.0;
//    double sum3 = 0.0;
//    for(int i=0; i<no; i++) {
//        for(int a=0; a<nv; a++) {
//            sum1 += hf.spinFock(i,a)*T1(i,a);           
//            for(int j=0; j<no; j++) {
//                for(int b=0; b<nv; b++) {
//                    sum2 += hf.spinTEI[i][j][a][b] * T2(i,j,a,b);
//                    sum3 += hf,spinTEI[i][j][a][b]*T1(i,a)*T(j,b);
//                }
//            }
//            e_cc += sum1;
//        }
//    }
//
//
//    return;
//}

//Delete allocated and used memory
HartreeFock::~HartreeFock()
{
    for(int i=0; i<2*norb; i++) {
        for(int j=0; j<2*norb; j++) {
            for(int k=0; k<2*norb; k++) {
                delete[] spinTEI[i][j][k];
            }
            delete[] spinTEI[i][j];
        }
        delete[] spinTEI[i];
    }
    delete[] spinTEI;

    for(int i=0; i<nmo; i++) {
        for(int j=0; j<nmo; j++) {
            for(int a=0; a<nmo; a++) {
                delete[] T2[i][j][a];
            }
            delete[] T2[i][j];
        }
        delete[] T2[i];
    }
    delete[] T2;

    for(int i=0; i<nmo; i++) {
        for(int j=0; j<nmo; j++) {
            for(int a=0; a<nmo; a++) {
                delete[] tau[i][j][a];
            }
            delete[] tau[i][j];
        }
        delete[] tau[i];
    }
    delete[] tau;

    for(int i=0; i<nmo; i++) {
        for(int j=0; j<nmo; j++) {
            for(int a=0; a<nmo; a++) {
                delete[] tau_t[i][j][a];
            }
            delete[] tau_t[i][j];
        }
        delete[] tau_t[i];
    }
    delete[] tau_t;
    
    for(int i=0; i<nmo; i++) {
        for(int j=0; j<nmo; j++) {
            for(int a=0; a<nmo; a++) {
                delete[] W_mnij[i][j][a];
            }
            delete[] W_mnij[i][j];
        }
        delete[] W_mnij[i];
    }
    delete[] W_mnij;
    
    for(int i=0; i<nmo; i++) {
        for(int j=0; j<nmo; j++) {
            for(int a=0; a<nmo; a++) {
                delete[] W_abef[i][j][a];
            }
            delete[] W_abef[i][j];
        }
        delete[] W_abef[i];
    }
    delete[] W_abef;


    for(int i=0; i<nmo; i++) {
        for(int j=0; j<nmo; j++) {
            for(int a=0; a<nmo; a++) {
                delete[] W_mbej[i][j][a];
            }
            delete[] W_mbej[i][j];
        }
        delete[] W_mbej[i];
    }
    delete[] W_mbej;

    for(int i=0; i<nmo; i++) {
        for(int j=0; j<nmo; j++) {
            for(int a=0; a<nmo; a++) {
                delete[] T2[i][j][a];
            }
            delete[] T2[i][j];
        }
        delete[] T2[i];
    }
    delete[] T2;

    for(int i=0; i<nmo; i++) {
        for(int j=0; j<nmo; j++) {
            for(int a=0; a<nmo; a++) {
                delete[] D_ijab[i][j][a];
            }
            delete[] D_ijab[i][j];
        }
        delete[] D_ijab[i];
    }
    delete[] D_ijab;

}

