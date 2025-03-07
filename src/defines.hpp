// defines.hpp
#pragma once

/* PARAMETRIZATION: Defines model of nuclear interactions */
#define USE_GM1
//#define USE_GM3
//#define USE_TM1

/* PARTICLES : Defines particles in the system */
//#define USE_NUCLEONS
#define USE_HYPERONS

/* MAGNETIC FIELD : Defines magnetic field in the system */
#define USE_MAGNETIC_FIELD
//#define USE_BI_MAGNETIC_FIELD
//#define USE_MODMAX_MAGNETIC_FIELD

/* LSV:  */
//#define LSV_TERMO_A
//#define LSV_TERMO_ISOLADO

/* BROYDEN: Defines Broyden tolerances and precision */
#define BROYDEN_EPS 1.0e-19
#define BROYDEN_TOLF 1.0e-12
#define BROYDEN_TOLX 1.0e-19
#define BROYDEN_STPMX 100.0
#define BROYDEN_MAXITS 40000
#define BROYDEN_NP 40
// FIX: NEED TO CHANGE BROYDEN'S FUNCTION IN ACCORDANCE WITH THE DEFINITIONS

// --------------------------------------------------------------------------------
// ERROR HANDLING

// Multiple definitions of parametrization
#if defined(GM1) + defined(GM3) + defined(TM1) > 1
#error "More than one parametrization is defined. Please define only one."
#endif

// Multiple definitions of particles
#if defined(NUCLEONS) + defined(HYPERONS) > 1
#error "More than one particle type is defined. Please define only one."
#endif

// Multiple definitions of magnetic field
#if defined(MAGNETIC_FIELD) + defined(BI_MAGNETIC_FIELD) + defined(MODMAX_MAGNETIC_FIELD) > 1
#error "More than one magnetic field type is defined. Please define only one."
#endif

// Multiple definitions of LSV
#if defined(LSV_TERMO_A) + defined(LSV_TERMO_ISOLADO) > 1
#error "More than one LSV type is defined. Please define only one."
#endif


// --------------------------------------------------------------------------------
// Define global variables
#define pi 3.14159265358979323846f
#define hc 197.3269804f      // hbar*c in MeV*fm
#define m 938.9187137765f    // Mass of nucleon (proton) in MeV
#define N_LEPTONS 2          // Number of leptons
#define ML {0.511f/m, 105.6583745f/m} // Masses of leptons (electron, muon) in MeV divided by nucleon mass
#define ms 400.0f/m          // Mass of sigma (scalar) meson in MeV divided by nucleon mass
#define mv 783.0f/m          // Mass of omega (vector) meson in MeV divided by nucleon mass
#define mr 770.0f/m          // Mass of rho (isovector) meson in MeV divided by nucleon mass
#define qe sqrt(4.0f * pi / 137.0f) // Electric charge in natural units
#define bc ml[0]**2/qe       // Critical magnetic field in natural units
#define alphaa 3.0f
#define betaa 1.0e-2f

// Parametrization
#ifdef USE_GM3
    #define n0 = 0.153f             // Saturation density in fm^-3
    #define gs sqrt(9.927f)/hc*m    // Couping constant for sigma
    #define gv sqrt(4.820f)/hc*m    // Couping constant for omega
    #define gr sqrt(4.791f)/hc*m    // Couping constant for rho
    #define rb 0.008659f            //
    #define rc 0.002421f            //
    #define rxi 0.0f                //
    #define xs 0.7f                 // 
    #define xv 0.783f               //
#elif defined(USE_GM1)
    #define n0 0.153f             // Saturation density in fm^-3
    #define gs sqrt(11.785f)/hc*m   // Couping constant for sigma
    #define gv sqrt(7.148f)/hc*m    // Couping constant for omega
    #define gr sqrt(4.41f)/hc*m     // Couping constant for rho
    #define rb 0.002948f            //
    #define rc 0.001071f            //
    #define rxi 0.0f                //
    #define xs 0.7f                 // 
    #define xv 0.783f               //
#elif defined(USE_TM1)  // FIX: NEED TO CHANGE TM1 PARAMETERS
    #define n0 = 0.145f         // Saturation density in fm^-3
    #define gs sqrt(11.785f)/hc*m   // Couping constant for sigma
    #define gv sqrt(7.148f)/hc*m    // Couping constant for omega
    #define gr sqrt(4.41f)/hc*m     // Couping constant for rho
    #define rb 0.002948f            //
    #define rc 0.001071f            //
    #define rxi 0.0f                //
    #define xs 0.7f                 // 
    #define xv 0.783f               //
#endif

// Particles
#if defined(USE_NUCLEONS)
    #define MB {939.56534623f/m, 1.0f} // Masses of Baryons in MeV divided by nucleon mass {neutron, proton}
    #define N_BARYONS 2
#elif defined(USE_HYPERONS)
    #define MB {939.56534623f/m, 1.0f, 1116.0f/m, 1193.0f/m, 1193.0f/m, 1193.0f/m, 1318.0f/m, 1318.0f/m} // Masses of Baryons in MeV divided by nucleon mass {neutron, proton, lambda0, sigma-, sigma0, sigma+, xi-, xi0}
    #define N_BARYONS 8
#endif

// Types definitions
typedef unsigned int uint;
typedef unsigned long ulong;