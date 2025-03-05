#include <setup.hpp>

void main_setup() {
    
    std::vector<double> x = {1.5, 1.5, 1.5}; // Initial guess

    std::cout << "Initial guess: ";
    for (double xi : x) std::cout << xi << " ";
    std::cout << std::endl;

    // Solve the system using Broyden's method
    broyden(x);

    // Display the solution
    std::cout << "Solution found: ";
    for (double xi : x) std::cout << xi << " ";
    std::cout << std::endl;
}