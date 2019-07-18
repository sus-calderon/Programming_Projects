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
    is >> natom;
    cout << "Number of atoms: " << natom << endl;

    //Now that I have natoms, I can create an natom X natom matrix
    S = new double* [natom];    //p2p2i
    for(int i=0; i<natom; i++)
        S[i] = new double[natom];   //We started with natom rows and now each row will have natom columns

    //Instead of doing a loop, I'm going to try to use the first two values in the line for my matrix
    //Now put values into space made
    //My file is like 28 lines so I need to loop 28 times to read each line
    is.clear();
    is.seekg(0, ios::beg);   // B/c I was at end of file, I need to clear flags I was at end, so I can 
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
        is >> m >> n >> S[m-1][n-1];
        start++;
    }

    is.close(); //Close input file
}

//Print out integrals we've read in
void HartreeFock::print_matrix()
{
    for(int i=0; i<natom; i++) {
        for(int j=0; j<natom; j++) {
            printf("%14.6f", S[i][j]);
        }
        printf("\n");
    }
}

//Delete allocated and used memory
HartreeFock::~HartreeFock()
{
    for(int i=0; i<natom; i++)
        delete[] S[i];
    delete[] S;
}


