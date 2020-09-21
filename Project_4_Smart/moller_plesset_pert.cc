// <moller_plesset_pert.cc>
// Written on August 29, 2019 by Susana Calderon
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cassert>
#include <cmath>
#include "hartreefock.h"
#include "molecule.h"
#include "moller_plesset_pert.h"
//Include Eigen package for easy diagonalization
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;  //Eigen::RowMaj    or makes it so that Eigen stores matrices by defaut in row-major order vs its actual default column-maj    or order
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

MP2::MP2(const char *filename)
{
    //Use this space to allocate data with Eigen package
}

//Transform the 2-electron Integrals to the MO Basis
int MP2::transform_AO_2_MO(HartreeFock& hf, MP2& mp2)
{

}

