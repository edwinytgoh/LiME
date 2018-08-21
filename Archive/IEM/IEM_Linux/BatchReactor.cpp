#include "cantera/thermo.h"
#include "cantera/zerodim.h"
#include "cantera/IdealGasMix.h"
#include <cstdlib.h>
#include <cstdio>
using namespace Cantera; 
class BatchReactor { 
    public: 
        BatchReactor(double T, double P, string X) {

        }
    private:
        ReactorNet *reactorNet;
        ConstPressureReactor *reactor; 
        IdealGasMix *gas;


};