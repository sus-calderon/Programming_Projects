//So this will be my .cc file that contains the definitions of my variables and functions
//I define stuff here.
//Written by Susana Calderon on June 17, 2019
// <hessian.cc>
//
#include "hessian.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cassert>  //Aborts program if user-specified condition NOT true
#include <cmath>    //Let's me use sqrt function

//For my Hessian function, I define how to obtain info from an input file
//This is also how I allocate memory
Hessian::Hessian(const char *filename)
{
    //Open filename here
    std::ifstream hessian(filename);     //"hessian" name of input file in code
                                    //filename generic variable name for file taken in as input
    assert(hessian.good());      //assert will tell me if a file was opened or not
    hessian >> natom;        //First thing it reads is the number of atoms in the molecule
    
    //Allocate space using next lines to be read
    //I need to keep reading in from my hessian file, but the number of lines of derivatives does not equal the number of atoms
    //double **H = new double* [natom*3];     
    H = new double* [natom*3];
    for(int i=0; i < natom*3; i++)   //Each atom has 3 coordinates thus natom*3
                                    //For water we have 3*3=9 so up to i=9
        H[i] = new double[natom*3]; //so my hessian has to be 3N by 3N so create a matrix that is 3N (will be square)

    //Now I need to fill up my matrix with values from the input hessian to be read
    for(int i=0; i < natom*3; i++) {
        for(int j=0; j < natom; j++) { // "i" is rows while "j" is column I believe
           // fscanf (hessian, "%lf %lf %lf", &H[i][3*j], &H[i][3*j+1], &H[i][3*j+2]);
           hessian >> H[i][3*j] >> H[i][3*j+1] >> H[i][3*j+2];
        }
    } 

    hessian.close();     //We then close input file when done.
}

//Now to print out my sqaure Hessian Matrix
void Hessian::print_H()
{
    for(int i=0; i<natom*3; i++) 
           printf("%20.12f %20.12f %20.12f %20.12f %20.12f %20.12f %20.12f %20.12f %20.12f\n", H[i][0], H[i][1], H[i][2], H[i][3], H[i][4], H[i][5], H[i][6], H[i][7], H[i][8]);
}

//And now for my deconstructor that deletes allocated memory.
Hessian::~Hessian()
{
    for(int i=0; i<natom*3; i++)
        delete[] H[i];
    delete[] H;
}
//Because I allocated memory at the start (manually) I must also deallocate/delete memory manually too


