// Ignition delay calculation with OpenMP. This example shows how to use OpenMP
// to run multiple reactor network calculations in parallel by using separate
// Cantera objects for each thread.

#include "cantera/zerodim.h"
#include "cantera/IdealGasMix.h"

#include <omp.h>

using namespace Cantera;

void arange(double start, double stop, double step, std::vector<double>& var)
{
    int N = (int) (stop - start)/step; 
    var.reserve(N);
    for (int i = 0; i < N; i++) {
        double current_entry = start + i*step;
        var.push_back(current_entry);
    }
}

void run()
{
    // The number of threads can be set by setting the OMP_NUM_THREADS
    // environment variable before running the code.
    int nThreads = omp_get_max_threads();
    writelog("Running on {} threads\n\n", nThreads);

    // Containers for Cantera objects to be used in different. Each thread needs
    // to have its own set of linked Cantera objects. Multiple threads accessing
    // the same objects at the same time will cause errors.
    std::vector<std::unique_ptr<IdealGasMix>> gases;
    std::vector<std::unique_ptr<IdealGasConstPressureReactor>> reactors;
    std::vector<std::unique_ptr<ReactorNet>> nets;

    // Create and link the Cantera objects for each thread. This step should be
    // done in serial
    for (int i = 0; i < nThreads; i++) {
        gases.emplace_back(new IdealGasMix("gri30.xml", "gri30"));
        reactors.emplace_back(new IdealGasConstPressureReactor());
        nets.emplace_back(new ReactorNet());
        reactors.back()->insert(*gases.back());
        nets.back()->addReactor(*reactors.back());
    }

    // Points at which to compute ignition delay time
    size_t nPoints = 50;
    vector_fp T0(nPoints);
    vector_fp ignition_time(nPoints);
    for (size_t i = 0; i < nPoints; i++) {
        T0[i] = 1000 + 500 * ((float) i) / ((float) nPoints);
    }

    // Calculate the ignition delay at each initial temperature using multiple
    // threads.
    //
    // Note on 'schedule(static, 1)':
    // This option causes points [0, nThreads, 2*nThreads, ...] to be handled by
    // the same thread, rather than the default behavior of one thread handling
    // points [0 ... nPoints/nThreads]. This helps balance the workload for each
    // thread in cases where the workload is biased, e.g. calculations for low
    // T0 take longer than calculations for high T0.
    #pragma omp parallel for schedule(static, 1)
    for (size_t i = 0; i < nPoints; i++) {
        // Get the Cantera objects that were initialized for this thread
        size_t j = omp_get_thread_num();
        IdealGasMix& gas = *gases[j];
        Reactor& reactor = *reactors[j];
        ReactorNet& net = *nets[j];

        // Set up the problem
        gas.setState_TPX(T0[i], OneAtm, "CH4:0.5, O2:1.0, N2:3.76");
        reactor.syncState();
        net.setInitialTime(0.0);

        // t = np.arange(0,estimatedIgnitionDelayTime,0.1*ms)
        
        // for i in range(len(t)):
        //     tnow = t[i]
        //     reactorNetwork.advance(tnow)
        //     timeHistory.loc[tnow] = reactorNetwork.get_state()
        // #print(timeHistory['OH'])

        // tau_ig = timeHistory['OH'].idxmax()

        // t1 = time.time()
        // runtime = t1-t0
        auto t = std::vector<double>();
        double end_time = 1.0; 
        double time_step = 0.01 * 1e-3;
        arange(0, end_time, time_step, t);

        for (auto t_now : t) {
            net.advance(t_now);
            
        }

    }

    // Print the computed ignition delays
    writelog("  T (K)    t_ig (s)\n");
    writelog("--------  ----------\n");
    for (size_t i = 0; i < nPoints; i++) {
        writelog("{: 8.1f}  {: 10.3e}\n", T0[i], ignition_time[i]);
    }
}

int main()
{
    try {
        run();
        appdelete();
        return 0;
    } catch (CanteraError& err) {
        // handle exceptions thrown by Cantera
        std::cout << err.what() << std::endl;
        std::cout << " terminating... " << std::endl;
        appdelete();
        return 1;
    }
}
