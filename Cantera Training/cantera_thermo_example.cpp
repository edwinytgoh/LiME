#include "cantera/thermo.h"
#include <iostream>
#include <cstdlib>

// using namespace Cantera

void thermo_demo(const std::string& file, const std::string& phase, double phi) {
    std::shared_ptr<Cantera::ThermoPhase> gas(Cantera::newPhase(file, phase)); 
    // std::stringstream X; // composition 
    // X << "CH4: " << phi << ", O2:" << 
    gas -> Cantera::ThermoPhase::setState_TPX(1500.0, 2.0*Cantera::OneAtm, "CH4:1, O2:2"); 
    // print T, P, rho
    std::cout << "------------------------------\n";
    std::cout << "Starting thermo_demo..." << std::endl;
    std::cout << "T = " << gas -> temperature() << " K" << std::endl; 
    std::cout << "P = " << gas -> pressure() << " Pa" << std::endl;
    std::cout << "rho = " << gas -> Cantera::ThermoPhase::density() << " kg/m3" << std::endl;

    // specific (per unit mass) thermodynamic properties
    std::cout << "Enthalpy per unit mass = " << gas -> enthalpy_mass() << " kJ/kg-K" << std::endl; 
    std::cout << "Entropy per unit mass = " << gas -> entropy_mass() << " kJ/kg-K" << std::endl;

    // chemical potentials of the species 
    int numSpecies = gas -> Cantera::ThermoPhase::nSpecies(); 
    std::cout << "Number of Species = " << numSpecies << std::endl;
    // std::unique_ptr<Cantera::vector_fp> mu(new Cantera::vector_fp(numSpecies)); 
    Cantera::vector_fp mu(numSpecies);
    gas -> getChemPotentials(&(mu[0])); 
    for (int n = 0; n < numSpecies; n++) {
        std::cout << "mu_" << gas -> speciesName(n) << " " << mu[n] << std::endl;
    }
    std::cout << "===============================================\n"; 
    std::cout << "Calculating equilibrium...\n";
    gas -> Cantera::ThermoPhase::equilibrate("HP"); 
    Cantera::
    std::cout << "Equilibrium T = " << gas -> temperature() << " K" << std::endl; 
}

int main(int argc, char **argv)
{
    std::unique_ptr<Cantera::ThermoPhase> gas(Cantera::newPhase("gri30.xml", "gri30"));
    std::cout << "T_init = " << gas->temperature() << " K" << std::endl;
    // try {
        thermo_demo("gri30.xml", "gri30", argv[1]);
    // } catch (CanteraError& err) {
    //     std::cout << err.what() << std::endl;
    //     return 1;
    // }
    return 0;
}