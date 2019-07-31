// <HF_SCF.cc>
// Written on July 16, 2019 by Susana Calderon
// This is my main program tht will use the declared and defined variables/functions in hartreefock.h and .cc
//

#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <cmath>
#include <cassert>
#include "hartreefock.h"
#include "molecule.h"

using namespace std;

int main()
{
    Molecule mol("geom.dat", 0);
    HartreeFock hf("s.dat");

    //Read and print nuclear repulsion energy
    hf.read_enuc("enuc.dat");

    //Read and print one electron integrals
    hf.S = hf.read_overlap(hf, "s.dat");
    hf.print_matrix("Overlap Integral Matrix (s): \n", hf.S);

    hf.T = hf.read_kinetic(hf, "t.dat");
    hf.print_matrix("Kinetic Energy Integral Matrix (t): \n", hf.T);

    hf.V = hf.read_potential(hf, "v.dat");
    hf.print_matrix("Nuclear Attraction Integral Matrix (v): \n", hf.V);

    hf.core = hf.build_core(hf);
    hf.print_matrix("Core Hamiltonian Matrix (h): \n", hf.core);

    //Read and print two electron integral
    hf.TEI = hf.read_tei(hf, "eri.dat");

    //Build Orthogonalization Matrix
    hf.SOM = hf.build_orthog(hf);
    hf.print_matrix("Symmetric Orthogonalization Matrix (S^1/2): \n", hf.SOM); 

    //Build Initial Guess Density
    //Build Initial (guess) Fock Matrix
    hf.F_Guess = hf.build_fock_guess(hf);
    hf.print_matrix("Initial Fock Matrix: \n", hf.F_Guess); 
    //Build Initial MO coefficient
    hf.MO_coef = hf.build_MO_coef(hf);
    hf.print_matrix("Initial Coefficient Matrix: \n", hf.MO_coef);
    //Build Density Matrix
    //hf.D = hf.build_density(hf.MO_coef);
    mol.electron = mol.electron_count();
    hf.build_density(hf.MO_coef, mol.electron);
    //hf.print_matrix("Initial Density Matrix: \n", hf.D);


    return 0;
}

