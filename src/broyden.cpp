#include "broyden.hpp"
#include <functional>

// Global parameters
constexpr int NP = BROYDEN_NP;         // Maximum system size
constexpr int MAXITS = BROYDEN_MAXITS;      // Maximum number of iterations
constexpr double EPS = BROYDEN_EPS;   // Numerical precision
constexpr double TOLF = BROYDEN_TOLF;   // Tolerance for the function
constexpr double TOLX = BROYDEN_TOLX;       // Tolerance for the solution
constexpr double STPMX = BROYDEN_STPMX;     // Maximum step size

// Define global variables
std::vector<double> fvec(NP);

// Objective function to minimize: 0.5 * ||f(x)||^2
double fmin(const std::vector<double>& x, std::function<std::vector<double>(const std::vector<double>&)> funcv) {
    fvec = funcv(x);
    double sum = 0.0;
    for (double val : fvec) {
        sum += val * val;
    }
    return 0.5 * sum;
}

// Gaussian elimination for solving linear systems (A * p = b)
std::vector<double> gaussElimination(const std::vector<std::vector<double>> &A_in, 
                                     const std::vector<double> &b_in) {
    int n = A_in.size();
    std::vector<std::vector<double>> A = A_in;
    std::vector<double> b = b_in;
    std::vector<double> x(n, 0.0);

    // Forward elimination
    for (int i = 0; i < n; ++i) {
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(A[k][i]) > std::abs(A[maxRow][i])) {
                maxRow = k;
            }
        }
        std::swap(A[i], A[maxRow]);
        std::swap(b[i], b[maxRow]);

        for (int k = i + 1; k < n; ++k) {
            double factor = A[k][i] / A[i][i];
            for (int j = i; j < n; ++j) {
                A[k][j] -= factor * A[i][j];
            }
            b[k] -= factor * b[i];
        }
    }

    // Back substitution
    for (int i = n - 1; i >= 0; --i) {
        x[i] = b[i] / A[i][i];
        for (int k = i - 1; k >= 0; --k) {
            b[k] -= A[k][i] * x[i];
        }
    }
    return x;
}

// Broyden’s method for solving f(x) = 0
std::vector<double> broyden(std::vector<double>& x, std::function<std::vector<double>(const std::vector<double>&)> funcv) {
    int n = x.size();

    std::vector<std::vector<double>> B(n, std::vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        B[i][i] = 1.0;  // Initial Jacobian approximation as identity matrix
    }

    std::vector<double> x_old(n), f_old(n), dx(n), df(n);
    fvec = funcv(x);

    for (int iter = 0; iter < MAXITS; ++iter) {
        f_old = fvec;
        x_old = x;

        // Solve B * dx = -fvec using Gaussian elimination
        for (int i = 0; i < n; ++i) {
            fvec[i] = -fvec[i];
        }
        dx = gaussElimination(B, fvec);

        // Update x
        for (int i = 0; i < n; ++i) {
            x[i] = x_old[i] + dx[i];
        }

        // Compute new function values
        fvec = funcv(x);
        if (std::sqrt(fmin(x, funcv)) < TOLF) {
            return fvec; // Converged successfully
        }

        // Compute df = f_new - f_old
        for (int i = 0; i < n; ++i) {
            df[i] = fvec[i] - f_old[i];
        }

        // Rank-1 update: B_new = B + ((df - B * dx) * dx^T) / (dx^T * dx)
        std::vector<double> Bdx(n);
        for (int i = 0; i < n; ++i) {
            Bdx[i] = 0.0;
            for (int j = 0; j < n; ++j) {
                Bdx[i] += B[i][j] * dx[j];
            }
        }

        double dx_norm = 0.0;
        for (double d : dx) {
            dx_norm += d * d;
        }
        if (dx_norm < EPS) {
            std::cerr << "Error: Division by zero."; // Avoid division by zero
        }

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                B[i][j] += (df[i] - Bdx[i]) * dx[j] / dx_norm;
            }
        }
    }

    std::cerr << "Warning: Broyden's method did not converge after " << MAXITS << " iterations.\n";
}
