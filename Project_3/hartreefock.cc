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

    //Now that I have natoms, I can create an natom X natom matrix
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


//Print out integrals we've read in
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
}


