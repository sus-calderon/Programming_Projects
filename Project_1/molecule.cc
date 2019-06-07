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
Molecule::Molecule(const char *filename, int q)
{
    charge = q;
    //Open filename here
    std::ifstream is(filename);     //"is" name of input file in code
                                    //filename generic variable name for file taken in as input
    assert(is.good());      //assert will tell me if a file was opened or not
    is >> natom;        //First thing it reads is the number of atoms in the molecule
    
    //Allocate space using next lines to be read

