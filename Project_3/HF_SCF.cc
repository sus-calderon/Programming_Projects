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

    //Read input integrals
    HartreeFock s("s.dat");
    cout << endl;
    cout << "Overlap Integral Matrix (s): \n";
    s.print_matrix();
    cout << endl;

    HartreeFock v("v.dat");
    cout << endl;
    cout << "Nuclear Attraction Integral Matrix (v): \n";
    v.print_matrix();
    cout << endl;

    HartreeFock t("t.dat");
    cout << endl;
    cout << "Kinetic Energy Integral Matrix (t): \n";
    t.print_matrix();
    cout << endl;

    return 0;
}

