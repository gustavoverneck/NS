#include "setup.hpp"

void main_setup() {
    NS neutron_star(0.1f, 1.50f, 0.01f, 1.0e15f);
    neutron_star.initialize();
    for (uint i = 0; i < neutron_star.npoints; ++i) {
        neutron_star.solve();
    }
}