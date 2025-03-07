#pragma once

#include "broyden.hpp"
#include "defines.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <array>

using namespace std;

class NS {
    public:
        // Constructors
        // Magnetized
        #if defined(USE_MAGNETIC_FIELD) && !defined(LSV_TERMO_A) && !defined(LSV_TERMO_ISOLADO)
        NS(const float mu_n_inf, const float mu_n_sup, const float mu_n_step, double B_G);
    #elif defined(USE_MAGNETIC_FIELD) && (defined(LSV_TERMO_A) || defined(LSV_TERMO_ISOLADO))
        NS(const float mu_n_inf, const float mu_n_sup, const float mu_n_step, double B_G, double xi);
    #elif (defined(USE_BI_MAGNETIC_FIELD) || defined(USE_MODMAX_MAGNETIC_FIELD)) && !defined(LSV_TERMO_A) && !defined(LSV_TERMO_ISOLADO)
        NS(const float mu_n_inf, const float mu_n_sup, const float mu_n_step, double B_G, double g);
    #elif (defined(USE_BI_MAGNETIC_FIELD) || defined(USE_MODMAX_MAGNETIC_FIELD)) && (defined(LSV_TERMO_A) || defined(LSV_TERMO_ISOLADO))
        NS(const float mu_n_inf, const float mu_n_sup, const float mu_n_step, double B_G, double g, double xi);
    #else
        NS(const float mu_n_inf, const float mu_n_sup, const float mu_n_step);
    #endif

        // Methods
        void initialize(); // Initialize all variables
        void solve();     // Solve the system of equations
        void output(string filename);    // Output the results

        // Variables
        uint npoints;   // Number of points in the grid

    private:
        // variables
        double mu_n_inf = 0.0f; // Neutron chemical potential lower limit
        double mu_n_sup = 0.0f; // Neutron chemical potential upper limit
        double mu_n_step = 0.0f; // Neutron chemical potential step size
        double B_G = 0.0f; // Magnetic field strength in Gauss
        double B_T = 0.0f; // Magnetic field in Tesla
        double B_surf = 1.0e11f; // Surface magnetic field in Gauss
        array<double, N_LEPTONS> ml = ML; // Leptons masses
        array<double, N_BARYONS> mb = MB; // Baryons masses
        double nuclear_magneton = 0.5f * qe / ml[0]; // Nuclear magneton in natural units
        double mu_n = 0.0f; // Neutron chemical potential
        double mu_e = 0.0f; // Electron chemical potential
        vector<double> fvec{4, 0.0f}; // Vector of nonlinear equations
        vector<double> x{4, 0.0f}; // Vector of variables {mu_e, sigma0, omega0, rho0}

        //File names
        string output_filename = "eos.csv";

        // Functions
        vector<double> funcv(const vector<double>& x); // Evaluate the system of nonlinear equations
        vector<double> map_x(const vector<double>& x); // Map x to the new variables
        float energy_density(const vector<double>& x); // Calculate the energy density
        float pressure(const vector<double>& x); // Calculate the pressure
        vector<double> eos = {0.0f, 0.0f}; // Equation of state {pressure, energy density}
    };