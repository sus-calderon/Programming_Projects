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
    hf.enuc = hf.read_enuc("enuc.dat");

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
    //hf.print_vector("TEI array: \n", hf.TEI);

    //Build Orthogonalization Matrix
    hf.SOM = hf.build_orthog(hf);
    hf.print_matrix("Symmetric Orthogonalization Matrix (S^1/2): \n", hf.SOM); 

    //Build Initial Guess Density
    //Build Initial (guess) Fock Matrix
    hf.F_Guess = hf.build_fock_guess(hf);
    hf.print_matrix("Initial Fock Matrix (F'): \n", hf.F_Guess); 
    //Build Initial MO coefficient
    hf.MO_coef = hf.build_MO_coef(hf);
    hf.print_matrix("Initial Coefficient Matrix (C): \n", hf.MO_coef);
    //Build Density Matrix
    hf.D = hf.build_density(hf, mol.electron_count());
    hf.print_matrix("Initial Density Matrix (D): \n", hf.D);

    //Compute the Initial SCF energy
    hf.SCF = hf.compute_SCF(hf);
    printf("The initial SCF electronic energy is %12.12f Hartrees.\n", hf.SCF);
    printf("The total energy (sum of the SCF electronic energy and nuclear repulsion energy) is %12.12f Hartrees.\n", (hf.SCF + hf.enuc)); 

    //Build the new Fock matrix (F) for the SCF proecure
    hf.F = hf.compute_Fock(hf);
    hf.print_matrix("Fock Matrix (F): \n", hf.F);

    return 0;
}

