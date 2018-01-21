// Command to compile and run: rm -rf batch_reactor.ex; g++ -g -o batch_reactor.ex batch_reactor.cpp $(pkg-config --cflags --libs cantera); ./batch_reactor.ex
#include "cantera/zerodim.h"
#include "cantera/thermo.h"
#include "cantera/IdealGasMix.h"
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <omp.h> 

int main()
{
    const std::string file = "gri30.xml";
    const std::string phase = "gri30";
    // std::shared_ptr<Cantera::ThermoPhase> gas(Cantera::newPhase(file, phase)); // NOTE: Cannot put ThermoPhase into a reactor. Whatever gets inserted into a reactor has to be a subclass of both Kinetics and ThermoPhase
    std::shared_ptr<Cantera::IdealGasMix> gas(new Cantera::IdealGasMix(file));

    // Cantera::ThermoPhase gas = *Cantera::newPhase(file, phase); // NOTE: This does not work because the TP.report function throws an error. This is a cantera bug. Also, the ThermoPhase() constructor should not be explicitly called.
    std::cout << "Creating Cantera IdealGasMix Object at address: " << gas.get() << "\n"
              << gas.get()->report() << std::endl;
    int T = 1000;
    int P = 25 * Cantera::OneAtm;
    double phi = 0.65;
    double dt = 0.0001, end_time = 10.; // seconds
    int num_points = (int) std::ceil(end_time / dt); 
    double* t = (double*) std::malloc(num_points * sizeof(double)); 
    const std::string outputFile = "br_output_withLimitedMaxStep.txt"; 
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < num_points; i++) {
        t[i] = i * dt; 
    }

    std::string X = "CH4:" + std::to_string(phi) + ",O2:" + std::to_string(2) + ",N2:" + std::to_string(2. * 0.79/0.21); // composition variable
    
    std::printf("X = %s\n", X.c_str()); 

    gas -> setState_TPX(T, P, std::string(X));
    // std::cout << gas -> report();
    
    Cantera::ConstPressureReactor br; 

    br.insert(*gas); // br needs to take the VALUE of gas, so use * to dereference. 
    std::cout << br.contents().report() << std::endl; 

    Cantera::ReactorNet reactor_network; 
    reactor_network.addReactor(br); 
    reactor_network.setMaxTimeStep(dt);
    std::ofstream outputFileStream; 
    outputFileStream.open(outputFile);
    for (int i = 0; i < num_points; i++) {
        std::cout << t[i] << std::endl;
        reactor_network.advance(t[i]); 
        outputFileStream << reactor_network.time() << " " << br.temperature() << " " << br.pressure() << " " << br.volume() << " " << br.density() << std::endl;
    }
    // reactor_network.advance(3600.0); 
    std::cout << "State of the reactor at t = " << std::to_string(reactor_network.time()) << " seconds " << std::endl;
    std::cout << reactor_network.reactor(0).contents().report();
    return 0;
}   