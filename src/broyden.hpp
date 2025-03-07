#pragma once

// Required includes for mathematical operations and vectors
#include "defines.hpp" // For global parameters
#include <iostream>  // For std::cout, std::cerr
#include <vector>    // For std::vector
#include <cmath>     // For std::abs, std::sqrt, etc.
#include <algorithm> // For std::max
#include <functional> // For std::function


// Global variables (extern to avoid multiple definition errors)
extern int nn;                    
extern std::vector<double> fvec;  

// Function prototypes

// Evaluates the system of nonlinear equations
std::vector<double> funcv(const std::vector<double>& x);

// Computes the objective function: 0.5 * ||f(x)||^2
double fmin(const std::vector<double>& x);

// Solves a linear system A * p = b using Gaussian elimination
std::vector<double> gaussElimination(const std::vector<std::vector<double>>& A_in,
                                     const std::vector<double>& b_in);

// Implements Broyden's method to solve f(x)=0
std::vector<double> broyden(std::vector<double>& x, std::function<std::vector<double>(const std::vector<double>&)> funcv);  