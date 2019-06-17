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
    zvals = new int[natom];     
    geom = new double* [natom]; //Making a pointer-to-pointer-to-int
                                // So a pointer points to a pointer that will point to my final #
    for(int i=0; i<natom; i++)
        geom[i] = new double[3];   // natom # of rows and then we set up coordinates next for each atom
        
    //Now we put things into the space made
    for(unsigned int i=0; i<natom; i++)
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

//This finds the unit vectors between atoms and their respective coordinates
// 0=x, 1=y, 2=z
double Molecule::unit(int cart, int a, int b)
{
    return -(geom[a][cart] -geom[b][cart])/bond(a,b);
    //This returns back the unit vector; it divides by the lenght between those atoms
}

//This function calculates the angles between atoms
double Molecule::angle(int i, int j, int k)
{
    return acos( unit(0, j, i)*unit(0, j, k)
                + unit(1, j, i)*unit(1, j, k)
                + unit(2, j, i)*unit(2, j, k) );
}

//This calculates the out of plane angle / dihedral angle between atoms
double Molecule::oop(int i, int j, int k, int l)
{
    double ejkl_x = ( unit(1,k,j)*unit(2,k,l) - unit(2,k,j)*unit(1,k,l) );
    double ejkl_y = ( unit(2,k,j)*unit(0,k,l) - unit(0,k,j)*unit(2,k,l) );
    double ejkl_z = ( unit(0,k,j)*unit(1,k,l) - unit(1,k,j)*unit(0,k,l) );

    double exx = ejkl_x * unit(0,k,i);
    double eyy = ejkl_y * unit(1,k,i);
    double ezz = ejkl_z * unit(2,k,i);

    double theta = (exx+eyy+ezz)/sin(angle(j,k,l));

    if(theta < -1.0) theta = asin(-1.0);
    else if(theta > 1.0) theta = asin(1.0);
    else theta = asin(theta);

    return theta;
}

//This calculates the torsion
//Computes the angle between planes a-b-c and b-c-d
// Computes the angle between planes a-b-c and b-c-d
double Molecule::torsion(int a, int b, int c, int d)
{
  double eabc_x = (unit(1,b,a)*unit(2,b,c) - unit(2,b,a)*unit(1,b,c));
  double eabc_y = (unit(2,b,a)*unit(0,b,c) - unit(0,b,a)*unit(2,b,c));
  double eabc_z = (unit(0,b,a)*unit(1,b,c) - unit(1,b,a)*unit(0,b,c));

  double ebcd_x = (unit(1,c,b)*unit(2,c,d) - unit(2,c,b)*unit(1,c,d));
  double ebcd_y = (unit(2,c,b)*unit(0,c,d) - unit(0,c,b)*unit(2,c,d));
  double ebcd_z = (unit(0,c,b)*unit(1,c,d) - unit(1,c,b)*unit(0,c,d));

  double exx = eabc_x * ebcd_x;
  double eyy = eabc_y * ebcd_y;
  double ezz = eabc_z * ebcd_z;

  double tau = (exx + eyy + ezz)/(sin(angle(a,b,c)) * sin(angle(b,c,d)));

  if(tau < -1.0) tau = acos(-1.0);
  else if(tau > 1.0) tau = acos(1.0);
  else tau = acos(tau);

  // Compute the sign of the torsion 
  double cross_x = eabc_y * ebcd_z - eabc_z * ebcd_y;
  double cross_y = eabc_z * ebcd_x - eabc_x * ebcd_z;
  double cross_z = eabc_x * ebcd_y - eabc_y * ebcd_x;
  double norm = cross_x*cross_x + cross_y*cross_y + cross_z*cross_z;
  cross_x /= norm;
  cross_y /= norm;
  cross_z /= norm;
  double sign = 1.0;
  double dot = cross_x*unit(0,b,c)+cross_y*unit(1,b,c)+cross_z*unit(2,b,c);
  if(dot < 0.0) sign = -1.0;

  return tau*sign;
}

//And now for my deconstructor that deletes allocated memory.
Molecule::~Molecule()
{
    delete[] zvals;
    for(int i=0; i<natom; i++)
        delete[] geom[i];
    delete[] geom;
}
//Because I allocated memory at the start (manually) I must also deallocate/delete memory manually too


