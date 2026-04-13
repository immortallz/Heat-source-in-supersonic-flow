#include "volterra.h"
#include <iostream>
#include <cmath>
#include <ctime>

int main() {
    // Initialize parameters
    VolterraParams params = get_default_volterra_params();

    // Print configuration
    std::cout << "Volterra solver for analytical solution" << std::endl;
    std::cout << "=======================================" << std::endl;
    std::cout << "Mach number: " << params.Mach << std::endl;
    std::cout << "Heat source intensity Q: " << params.Q << std::endl;
    std::cout << "Source position (r, x): (" << params.r_source << ", " << params.x_source << ")" << std::endl;
    std::cout << "Computation length L: " << params.L << std::endl;
    std::cout << "Grid points N: " << params.N << std::endl;
    std::cout << "=======================================" << std::endl;

    // Solve
    clock_t start = clock();
    VolterraResult result = solve_volterra(params);
    clock_t end = clock();

    // Save results
    save_results(result, params);

    // Print some results
    double h = params.L / params.N;
    std::cout << "Derivative at x=0: " << PI * (result.dphi1_dx[1] + result.dphiu1_dx[1]
                                               - result.dphi1_dx[0] - result.dphiu1_dx[0]) / h << std::endl;
    std::cout << "Execution time: " << static_cast<double>(end - start) / CLOCKS_PER_SEC << "s" << std::endl;

    std::cout << "Results saved to " << params.output_dir << "/" << std::endl;

    return 0;
}
