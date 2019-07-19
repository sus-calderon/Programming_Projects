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

using namespace std;

int main()
{
    //Read in enuc.dat file
    double enuc;
    ifstream nucl;
    nucl.open("enuc.dat");
    assert(nucl.good());
    nucl >> enuc;
    cout << "Nuclear Repulsion Energy: " << enuc << endl;
    cout << endl;

    HartreeFock hf("s.dat");

    //Read and print input integrals
    hf.read_oei(hf.S, "s.dat");
    hf.print_matrix("Overlap Integral Matrix (s): \n", hf.S);

    hf.read_oei(hf.T, "t.dat");
    hf.print_matrix("Kinetic Energy Integral Matrix (t): \n", hf.T);

    hf.read_oei(hf.V, "v.dat");
    hf.print_matrix("Potential Energy Integral Matrix (v): \n", hf.V);

    hf.build_core(hf.T, hf.V);
    hf.print_matrix("Core Hamiltonian Matrix (h): \n", hf.core);

    return 0;
}

