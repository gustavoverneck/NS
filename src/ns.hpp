#include "broyden.hpp"
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <array>

using namespace std;

class NS {
    public:
        #ifdef USE_MAGNETIC_FIELD
            NS(const int mu_n_inf, const int mu_n_sup, const int mu_n_step, double B_G); // Constructor
        #else
            NS(const int mu_b_inf, const int mu_b_sup, const int mu_b_step); // Constructor
        #endif

        // Methods
        void initialize(); // Initialize all variables
        void solve();     // Solve the system of equations
        void output(string filename);    // Output the results

    private:
        // variables
        double mu_n_inf = 0.0f; // Neutron chemical potential lower limit
        double mu_n_sup = 0.0f; // Neutron chemical potential upper limit
        double mu_n_step = 0.0f; // Neutron chemical potential step size
        uint npoints;   // Number of points in the grid
        double B_G = 0.0f; // Magnetic field strength in Gauss
        double B_T = 0.0f; // Magnetic field in Tesla
        double B_surf = 1.0e11f; // Surface magnetic field in Gauss
        array<double, N_LEPTONS> ml = ML; // Leptons masses
        array<double, N_BARYONS> mb = MB; // Baryons masses
        double nuclear_magneton = 0.5f * qe / ml[0]; // Nuclear magneton in natural units
        double mu_n = 0.0f; // Neutron chemical potential
        double mu_e = 0.0f; // Electron chemical potential

        //File names
        string output_filename = "eos.csv";

        // Functions
        vector<double> funcv(const vector<double>& x); // Evaluate the system of nonlinear equations
        vector<double> map_x(const vector<double>& x); // Map x to the new variables
    };