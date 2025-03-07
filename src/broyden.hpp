#pragma once

// Required includes for mathematical operations and vectors
#include "defines.hpp" // For global parameters
#include <iostream>  // For std::cout, std::cerr
#include <vector>    // For std::vector
#include <cmath>     // For std::abs, std::sqrt, etc.
#include <algorithm> // For std::max
#include <functional> // For std::function

using namespace std;

// Global variables (extern to avoid multiple definition errors)
extern int nn;                    
extern vector<double> fvec;  

// Function prototypes

// Evaluates the system of nonlinear equations
vector<double> funcv(const vector<double>& x);

// Computes the objective function: 0.5 * ||f(x)||^2
double fmin(const vector<double>& x);

// Solves a linear system A * p = b using Gaussian elimination
vector<double> gaussElimination(const vector<vector<double>>& A_in, const vector<double>& b_in);

// Implements Broyden's method to solve f(x)=0
vector<double> broyden(vector<double>& x, function<vector<double>(const vector<double>&)> funcv);  