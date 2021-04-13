// <MP2.cc>
// Written on August 13, 2019 by Susana Calderon
// Main Program: The second-order Moller-Plesset perturbation (MP2) energy
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
    HartreeFock hf("s.dat", mol.electron_count());

    //Read and print nuclear repulsion energy
    hf.enuc = hf.read_enuc("enuc.dat");

    //Read and print one electron integrals
    hf.read_overlap(hf, "s.dat");
    hf.read_kinetic(hf, "t.dat");
    hf.read_potential(hf, "v.dat");
    hf.build_core(hf);

    //Read and print two electron integral
    hf.read_tei(hf, "eri.dat");

    //Build Orthogonalization Matrix
    hf.build_orthog(hf);

    //Build Initial Guess Density with guess Fock Matrix
    hf.F = hf.core;
    hf.build_density(hf, mol.electron_count());
    hf.old_D = hf.D;

    //Compute the Initial SCF energy
    hf.compute_SCF(hf);
    hf.old_SCF = hf.SCF;

    //Build the new Density matrix and iterate Step 7-9 until convergence is reached
    double tol = 1e-12;
    hf.iter_max = 100;
    hf.iter = 0;
    double delta_E = 1.0;
    while(abs(delta_E) > tol && hf.iter < hf.iter_max) {
        hf.iter+=1;
        double rms_D = 0.0;

        //Build the new Fock matrix (F) for the SCF proecure
        hf.update_Fock(hf);
        //Build the new Density matrix (D)
        hf.build_density(hf, mol.electron_count());
        //Compute the new SCF energy
        hf.compute_SCF(hf);

        //Test the SCF electronic energy for convergence
        delta_E = hf.SCF - hf.old_SCF;
        hf.old_SCF = hf.SCF;

        //Test the root-mean-squared difference in Densities for convergence
        for(int i=0; i<hf.D.rows(); i++) {
            for(int j=0; j<hf.D.cols(); j++) {
                rms_D += (hf.D(i,j) - hf.old_D(i,j))*(hf.D(i,j) - hf.old_D(i,j));
            }
        }
        rms_D = sqrt(rms_D); 
        hf.old_D = hf.D;

        //Print out iterations and corresponding SCF and total electronic energy, along with Delta E and RMS_D
        if(hf.iter == 1) {
            //printf("%-12s %-20s %-20s %-20s %-20s \n", "Iter", "E(elec)", "E(tot)", "Delta(E)", "RMS(D)");
        }
        //printf("%04d %20.12f %20.12f %20.12f %20.12f \n", hf.iter, hf.SCF, hf.tot_E, delta_E, rms_D);

    }
    cout << "HF SCF Energy = " << hf.tot_E << endl;

    //Transform the 2e- integrals into the MO basis
    hf.smart_transform_AO_2_MO(hf);

    //Calculate the MP2 energy using the TEI in the MO basis
    hf.smart_MP2_calc(hf, mol.electron_count());
    cout << "MP2 Energy Correction = " << hf.Emp2 << endl;
    
    cout << "Total Energy = " << hf.tot_E + hf.Emp2 << endl;

    //Print out converged Fock matrix
    hf.print_matrix("The converged Fock matrix in the AO basis:\n", hf.F);

    //CCSD Functions
    hf.transform_to_spin(hf);

    hf.create_spinFock(hf, mol.electron_count());

    cout << endl;

    hf.build_cluster_amp(hf, mol.electron_count());

    hf.build_tau(hf);

    hf.ccF_intermediates(hf, mol.electron_count());

    hf.ccW_intermediates(hf, mol.electron_count());

    hf.update_t_ia(hf);
    
    hf.update_t_ijab(hf);

    return 0;
}
