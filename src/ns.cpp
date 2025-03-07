#include "ns.hpp"

// Constructor for the NS class
#if defined(USE_MAGNETIC_FIELD) && !defined(LSV_TERMO_A) && !defined(LSV_TERMO_ISOLADO)
    NS::NS(const float mu_n_inf, const float mu_n_sup, const float mu_n_step, double B_G)
        : mu_n_inf(mu_n_inf), mu_n_sup(mu_n_sup), mu_n_step(mu_n_step), B_G(B_G) {}
#elif defined(USE_MAGNETIC_FIELD) && (defined(LSV_TERMO_A) || defined(LSV_TERMO_ISOLADO))
    NS::NS(const float mu_n_inf, const float mu_n_sup, const float mu_n_step, double B_G, double xi)
        : mu_n_inf(mu_n_inf), mu_n_sup(mu_n_sup), mu_n_step(mu_n_step), B_G(B_G), xi(xi) {}
#elif (defined(USE_BI_MAGNETIC_FIELD) || defined(USE_MODMAX_MAGNETIC_FIELD)) && !defined(LSV_TERMO_A) && !defined(LSV_TERMO_ISOLADO)
    NS::NS(const float mu_n_inf, const float mu_n_sup, const float mu_n_step, double B_G, double g)
        : mu_n_inf(mu_n_inf), mu_n_sup(mu_n_sup), mu_n_step(mu_n_step), B_G(B_G), g(g) {}
#elif (defined(USE_BI_MAGNETIC_FIELD) || defined(USE_MODMAX_MAGNETIC_FIELD)) && (defined(LSV_TERMO_A) || defined(LSV_TERMO_ISOLADO))
    NS::NS(const float mu_n_inf, const float mu_n_sup, const float mu_n_step, double B_G, double g, double xi)
        : mu_n_inf(mu_n_inf), mu_n_sup(mu_n_sup), mu_n_step(mu_n_step), B_G(B_G), g(g), xi(xi) {}
#else
    NS::NS(const float mu_n_inf, const float mu_n_sup, const float mu_n_step)
        : mu_n_inf(mu_n_inf), mu_n_sup(mu_n_sup), mu_n_step(mu_n_step) {}
#endif


void NS::initialize() {
    // Initialize all variables
    const uint npoints = (mu_n_sup - mu_n_inf) / mu_n_step;
    B_T = B_G / 1.0e4f; // Convert Gauss to Tesla
    vector<double> fvec{4, 0.0f}; // Initialize fvec
    vector<double> x(npoints, 0.0f); // Initialize x {mu_e, sigma0, omega0, rho0}
}

vector<double> NS::funcv(const vector<double>& x) {
    vector<double> f = {0.0f, 0.0f, 0.0f, 0.0f};
    f[0] = x[0] - mu_e; // Electron chemical potential
    f[1] = x[1] - (x[0] * x[0] + x[2] * x[2] + x[3] * x[3]) / (2.0f * n0); // Sigma field
    f[2] = x[2] - (x[0] * x[0] + x[1] * x[1] + x[3] * x[3]) / (2.0f * n0); // Omega field
    f[3] = x[3] - (x[0] * x[0] + x[1] * x[2] + x[2] * x[3]) / (2.0f * n0); // Rho field
    return f;
}

void NS::solve() {
    auto funcv_wrapper = [this](const std::vector<double>& x) {
        return this->funcv(x);
    };

    std::vector<double> x = broyden(x, funcv_wrapper);
}


vector<double> NS::map_x(const vector<double>& x) {
    // Map x to the new variables
    return x;
}


float NS::energy_density(const vector<double>& x) {
    // Calculate the energy density
    float e = 0.0f;
    return e;
}

float NS::pressure(const vector<double>& x) {
    // Calculate the pressure
    float p = 0.0f;
    return p;
}

void NS::output(string filename) {

}