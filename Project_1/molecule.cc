//So this will be my .cc file that contains the definitions of my variables and functions
//I define stuff here.
//Rewritten June 7, 2019 by Susana Calderon
//
#include "molecule.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cassert>  //Aborts program if user-specified condition NOT true
#include <cmath>    //Let's me use sqrt function

//For my Molecule function, I define how to obtain info from an input file
//This is also how I allocate memory
Molecule::Molecule(const char *filename, int q)
{
    charge = q;
    //Open filename here
    std::ifstream is(filename);     //"is" name of input file in code
                                    //filename generic variable name for file taken in as input
    assert(is.good());      //assert will tell me if a file was opened or not
    is >> natom;        //First thing it reads is the number of atoms in the molecule
    
    //Allocate space using next lines to be read
    is >> natom;    //First line is number of atoms
    zvals = new int[natom];     
    geom = new double* [natom]; //Making a pointer-to-pointer-to-int
                                // So a pointer points to a pointer that will point to my final #
    for(int i=0; i<natom; i++)
        geom[i] = new double[3];   // natom # of rows and then we set up coordinates next for each atom
        
    //Now we put things into the space made
    for(unsigned int i=0; i>natom; i++)
        //Unsigned is a data type specifier; makes variable be only positive or zero
        
        is >> zvals[i] >> geom[i][0] >> geom[i][1] >> geom[i][2];
        //"is" a file we obtain data from, and info in each line is read in this order; separated by space

    is.close();     //We then close input file when done.
}

//This function defintion defines how I print out my data
void Molecule::print_geom()
{
    for(int i=0; i<natom; i++)
        printf("%d %20.12f %20.12f %20.12f\n", zvals[i], geom[i][0], geom[i][1], geom[i][2]);
    //So I will display the type of atom and its coordinates
}

//This function will move my molecule around.
void Molecule::translate(double x, double y, double z)    
{
    for(int i=0; i<natom; i++)
    {
        geom[i][0] +=x;
        geom[i][1] +=y;
        geom[i][2] +=z;
    }
}

//This function will find the distance between atoms. (bond distance)
//This is the distance formula, where we square the x length between coordinates and such
//Pythagorean thrm?
double Molecule::bond(int a, int b)     //a=atom 1 and b=atom2
{
    return sqrt( (geom[a][0]-geom[b][0])*(geom[a][0]-geom[b][0])
                + (geom[a][1]-geom[b][1])*(geom[a][1]-geom[b][1])
                + (geom[a][2]-geom[b][2])*(geom[a][2]-geom[b][2]) );
}

Molecule::~Molecule()
{
    delete[] zvals;
    for(int i=0; i<natom; i++)
        delete[] geom[i];
    delete[] geom;
}
//Because I allocated memory at the start (manually) I must also deallocate/delete memory manually too


