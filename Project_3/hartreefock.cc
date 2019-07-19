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

void HartreeFock::read_oei(double** oei_mat, const char *filename)
{
    //Open File here
    std::ifstream is(filename);
    assert(is.good());
   
    //Instead of doing a loop, I'm going to try to use the first two values in the line for my matrix
    //Now put values into space made
    //My file is like 28 lines so I need to loop 28 times to read each line
    int lines=0;
    std::string line;
    while (getline(is, line)) {
       lines++;
    }
    //So it does read in 28 lines for s.dat

    is.clear();
    is.seekg(0, ios::beg);
    int start=1;
    int m;
    int n;
    while( start<=lines ) {
        is >> m >> n >> oei_mat[m-1][n-1];
        oei_mat[n-1][m-1] = oei_mat[m-1][n-1];
        start++;
    }

    is.close(); //Close input file

    return;
}


//Print out integrals we've read in
void HartreeFock::print_matrix(std::string mat_string, double** matrix)
{
    cout << endl;
    cout << mat_string;
    for(int i=0; i<norb; i++) {
        for(int j=0; j<norb; j++) {
            printf("%14.6f", matrix[i][j]);
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


