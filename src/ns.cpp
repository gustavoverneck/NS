#include "ns.hpp"

#ifdef USE_MAGNETIC_FIELD
NS::NS(const int mu_n_inf, const int mu_n_sup, const int mu_n_step, double B_G)
    :mu_n_inf(mu_n_inf), mu_n_sup(mu_n_sup), mu_n_step(mu_n_step), B_G(B_G) {
}
#else
NS::NS(const int mu_b_inf, const int mu_b_sup, const int mu_b_step)
    :mu_n_inf(mu_n_inf), mu_n_sup(mu_n_sup), mu_n_step(mu_n_step) {
}
#endif


void NS::initialize() {
    // Initialize all variables
    const uint npoints = (mu_n_sup - mu_n_inf) / mu_n_step;
    B_T = B_G / 1.0e4f; // Convert Gauss to Tesla
    vector<double> fvec{4, 0.0f}; // Initialize fvec
    vector<double> x(npoints, 0.0f); // Initialize x {mu_e, sigma0, omega0, rho0}
}


void NS::solve() {
    for (int i = 0; i < npoints; ++i) {
        vector<double> x = broyden(x);
    }
}


vector<double> NS::funcv(const vector<double>& x) {
    // Evaluate the system of nonlinear equations
}


vector<double> NS::map_x(const vector<double>& x) {
    // Map x to the new variables
}


void NS::output(string filename) {}