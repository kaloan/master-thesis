#include <algorithm>
#include <array>
#include <bits/floatn-common.h>
#include <chrono>
// #include <cmath>
#include <filesystem>
#include <format>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
// #include <iterator>
#include <memory>
#include <sys/cdefs.h>
#include <vector>
using namespace std;

// #undef y1
constexpr double M_PI = 3.141592653589793115997963468544185161590576171875;

template <typename T, size_t n, size_t m>
using Array2DGeneric = array<array<T, m>, n>;

template <size_t n, size_t m>
// using Array2D = array<array<array<array<double, n4>, n3>, m>, n>;
using Array2D = Array2DGeneric<double, n, m>;

// humans
constexpr auto N = 5.38e6;
constexpr auto Gamma = 0.071;

// mosquito
constexpr auto M = 5.38e7;
constexpr auto mu = 0.033;
constexpr auto a = 0.241;

constexpr auto BetaHV = 0.5;
constexpr auto BetaVH = 0.1;
constexpr auto tau = 10.;

// multipliers in normalized model
const auto b = BetaVH * __builtin_exp(-mu * tau) * a * M / N;
constexpr auto c = BetaHV * a;

// maximum infected
constexpr auto Imax = .1;
// set to reduce the size of the domain where the viability kernel lies
constexpr auto xmax = .1;
constexpr auto ymax = .2;

// control
constexpr auto k = .6;
constexpr auto umax = .6;

// lambda - lipschitz constant for the system
const auto Lexact = 3. * (b + c) + Gamma + mu;
const auto L = 1.05 * Lexact;

// numerical space step
constexpr auto Nx = 101;
constexpr auto dx = xmax / (Nx - 1);
constexpr auto Ny = 101;
constexpr auto dy = ymax / (Ny - 1);

// x and y
constexpr auto padX = 6;
array<double, Nx + 2 * padX> x;
constexpr auto padY = 6;
array<double, Ny + 2 * padY> y;
const auto hx = make_unique<Array2D<x.size(), y.size()>>();
const auto hy = make_unique<Array2D<x.size(), y.size()>>();
const auto mx = make_unique<Array2D<x.size(), y.size()>>();
const auto my = make_unique<Array2D<x.size(), y.size()>>();
const auto distSgnd = make_unique<Array2D<x.size(), y.size()>>();

const auto Fxp = make_unique<Array2D<x.size(), y.size()>>();
const auto Fxm = make_unique<Array2D<x.size(), y.size()>>();
const auto Fyp = make_unique<Array2D<x.size(), y.size()>>();
const auto Fym = make_unique<Array2D<x.size(), y.size()>>();

const auto wx = make_unique<Array2D<x.size() + 5, y.size()>>();
const auto wy = make_unique<Array2D<x.size(), y.size() + 5>>();

const auto hamiltonian = make_unique<Array2D<x.size(), y.size()>>();

template <size_t n, size_t m>
__always_inline void weightedSum(const unique_ptr<Array2D<n, m>> &res,
                                 const unique_ptr<Array2D<n, m>> &mat1,
                                 const unique_ptr<Array2D<n, m>> &mat2,
                                 const double &w1, const double &w2) noexcept {
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < m; ++j) {
      (*res)[i][j] = w1 * (*mat1)[i][j] + w2 * (*mat2)[i][j];
    }
  }
}

template <size_t n, size_t m>
__always_inline void printSolutionOmega(const Array2D<n, m> &arr,
                                        ostream &stream = std::cout,
                                        const string &delim = ", ") noexcept {
  for (size_t i = 0; i < Nx; ++i) {
    for (size_t j = 0; j < Ny; ++j) {
      stream << fixed << setprecision(8);
      stream << x[i + padX] << delim << y[j + padY] << delim;
      stream << scientific;
      stream << arr[i + padX + 3][j + padY + 3] << "\n";
    }
  }
}

template <size_t n, size_t m>
__always_inline void
printSolutionExtended(const Array2D<n, m> &arr, ostream &stream = std::cout,
                      const string &delim = ", ") noexcept {
  for (size_t i = 0; i < x.size(); ++i) {
    for (size_t j = 0; j < y.size(); ++j) {
      stream << fixed << setprecision(8);
      stream << x[i] << delim << y[j] << delim;
      stream << scientific;
      stream << arr[i + 3][j + 3] << "\n";
    }
  }
}

template <size_t n, size_t m>
__always_inline void
matrixContinuation(const unique_ptr<Array2D<n, m>> &mat) noexcept {
  for (size_t i = 0; i < n - 3; ++i) {
    (*mat)[i + 3][0] = (*mat)[i + 3][1] = (*mat)[i + 3][2] = (*mat)[i + 3][3];
    (*mat)[i + 3][m - 1] = (*mat)[i + 3][m - 2] = (*mat)[i + 3][m - 3] =
        (*mat)[i + 3][m - 4];
  }

  for (size_t j = 0; j < m - 3; ++j) {
    (*mat)[0][j + 3] = (*mat)[1][j + 3] = (*mat)[2][j + 3] = (*mat)[3][j + 3];
    (*mat)[n - 1][j + 3] = (*mat)[n - 2][j + 3] = (*mat)[n - 3][j + 3] =
        (*mat)[n - 4][j + 3];
  }

  for (size_t i = 0; i < 3; ++i) {
    for (size_t j = 0; j < 3; ++j) {
      (*mat)[i][j] = (*mat)[3][3];

      (*mat)[n - 3 + i][j] = (*mat)[n - 4][3];

      (*mat)[i][m - 3 + j] = (*mat)[3][m - 4];
      (*mat)[n - 3 + i][m - 3 + j] = (*mat)[n - 4][m - 4];
    }
  }
}

__always_inline double weno(const double &v1, const double &v2,
                            const double &v3, const double &v4,
                            const double &v5) noexcept {
  // Smoothness of stencils
  const auto S1 = 13. / 12. * (v1 - 2. * v2 + v3) * (v1 - 2. * v2 + v3) +
                  1. / 4. * (v1 - 4. * v2 + 3. * v3) * (v1 - 4. * v2 + 3. * v3);
  const auto S2 = 13. / 12. * (v2 - 2. * v3 + v4) * (v2 - 2. * v3 + v4) +
                  1. / 4. * (v2 - v4) * (v2 - v4);
  const auto S3 = 13. / 12. * (v3 - 2. * v4 + v5) * (v3 - 2. * v4 + v5) +
                  1. / 4. * (3. * v3 - 4. * v4 + v5) * (3. * v3 - 4. * v4 + v5);

  // Find the maximum of vi's
  const auto maxv2 = max({v1 * v1, v2 * v2, v3 * v3, v4 * v4, v5 * v5});
  const auto eps1 = 1e-6 * maxv2 + 1e-99;

  // Weights
  const auto a1 = .1 / ((S1 + eps1) * (S1 + eps1));
  const auto a2 = .6 / ((S2 + eps1) * (S2 + eps1));
  const auto a3 = .3 / ((S3 + eps1) * (S3 + eps1));
  const auto w1 = a1 / (a1 + a2 + a3);
  const auto w2 = a2 / (a1 + a2 + a3);
  const auto w3 = a3 / (a1 + a2 + a3);
  const auto f1x = v1 / 3. - 7. * v2 / 6. + 11. * v3 / 6.;
  const auto f2x = -v2 / 6. + 5. * v3 / 6. + v4 / 3.;
  const auto f3x = v3 / 3. + 5. * v4 / 6. - v5 / 6.;

  return w1 * f1x + w2 * f2x + w3 * f3x;
}

template <size_t n, size_t m>
__always_inline void
approximate_derivatives(const unique_ptr<Array2D<n, m>> &mat) noexcept {

  // Compute step forward matrices
  for (size_t i = 0; i < n - 1; ++i) {
    for (size_t j = 0; j < m - 6; ++j) {
      (*wx)[i][j] = ((*mat)[i + 1][j + 3] - (*mat)[i][j + 3]) / dx;
    }
  }
  for (size_t i = 0; i < n - 6; ++i) {
    for (size_t j = 0; j < m - 1; ++j) {
      (*wy)[i][j] = ((*mat)[i + 3][j + 1] - (*mat)[i + 3][j]) / dy;
    }
  }

  // Compute forward and backward approximations of derivatives
  for (size_t i = 0; i < n - 6; ++i) {
    for (size_t j = 0; j < m - 6; ++j) {
      (*Fxm)[i][j] = weno((*wx)[i][j], (*wx)[i + 1][j], (*wx)[i + 2][j],
                          (*wx)[i + 3][j], (*wx)[i + 4][j]);
      (*Fxp)[i][j] = weno((*wx)[i + 5][j], (*wx)[i + 4][j], (*wx)[i + 3][j],
                          (*wx)[i + 2][j], (*wx)[i + 1][j]);

      (*Fym)[i][j] = weno((*wy)[i][j], (*wy)[i][j + 1], (*wy)[i][j + 2],
                          (*wy)[i][j + 3], (*wy)[i][j + 4]);
      (*Fyp)[i][j] = weno((*wy)[i][j + 5], (*wy)[i][j + 4], (*wy)[i][j + 3],
                          (*wy)[i][j + 2], (*wy)[i][j + 1]);
    }
  }
}

__always_inline double hamiltonianInDirection(const double &fp,
                                              const double &fm, const double &h,
                                              const double &a) noexcept {
  return h * (fp + fm) / 2. - a * (fp - fm) / 2.;
}

__always_inline void calcHamiltonian() noexcept {

  for (size_t i = 0; i < x.size(); ++i) {
    for (size_t j = 0; j < y.size(); ++j) {
      auto hamilt = 0.;

      // Find whether we take the maximum
      const auto chooseUMax =
          ((*mx)[i][j] * ((*Fxp)[i][j] + (*Fxm)[i][j]) / 2. +
           (*my)[i][j] * ((*Fyp)[i][j] + (*Fym)[i][j]) / 2.) > 1e-99;

      // Partial derivatives respecting the directions
      const auto hxa = (*hx)[i][j] + chooseUMax * (*mx)[i][j];
      const auto hya = (*hy)[i][j] + chooseUMax * (*my)[i][j];

      const auto ax = max(abs((*hx)[i][j]), abs(hxa));
      const auto ay = max(abs((*hy)[i][j]), abs(hya));

      // Lax-Friedrichs step
      hamilt += hamiltonianInDirection((*Fxp)[i][j], (*Fxm)[i][j], hxa, ax);
      hamilt += hamiltonianInDirection((*Fyp)[i][j], (*Fym)[i][j], hya, ay);

      (*hamiltonian)[i][j] = hamilt;
    }
  }
}

template <size_t n, size_t m>
__always_inline void performTimeStep(const unique_ptr<Array2D<n, m>> &tp1Matrix,
                                     const unique_ptr<Array2D<n, m>> &tMatrix,
                                     const double &dt) noexcept {
  approximate_derivatives(tMatrix);
  calcHamiltonian();

  for (size_t i = 0; i < n - 6; ++i) {
    for (size_t j = 0; j < m - 6; ++j) {
      (*tp1Matrix)[i + 3][j + 3] = max((1 - L * dt) * (*tMatrix)[i + 3][j + 3] -
                                           dt * (*hamiltonian)[i][j],
                                       (*distSgnd)[i][j]);
    }
  }

  matrixContinuation(tp1Matrix);
}

__always_inline double signed_distance(const double &x, const double &y) {
  if (y >= 0) {
    if (x >= 0) {
      if (x <= Imax)
        return -abs(x - Imax);
      else
        return abs(x - Imax);
    } else
      return abs(x);
  } else {
    if (x >= 0) {
      if (x <= Imax)
        return abs(y);
      else
        return __builtin_sqrt((x - Imax) * (x - Imax) + y * y);
    } else
      return __builtin_sqrt(x * x + y * y);
  }
}

__always_inline void initialize() {

  // Fill gridpoint coordinates
  generate(x.begin(), x.end(),
           [n = 0]() mutable { return (-padX + n++) * dx; });
  x[padX] = 0;
  x[padX + Nx - 1] = xmax;
  generate(y.begin(), y.end(),
           [n = 0]() mutable { return (-padY + n++) * dy; });
  y[padY] = 0;
  y[padY + Ny - 1] = ymax;

  for (size_t i = 0; i < x.size(); ++i) {
    for (size_t j = 0; j < y.size(); ++j) {
      // Compute the partial derivative constant terms
      (*hx)[i][j] = Gamma * x[i] - (1. - x[i]) * b * y[j];
      (*hy)[i][j] = mu * y[j] - (1. - y[j]) * c * x[i];

      // Compute the partial derivative terms inside maximums
      (*mx)[i][j] = k * umax * (1. - x[i]) * b * y[j];
      (*my)[i][j] = k * (1. - y[j]) * c * umax * x[i];

      // Compute the signed distance funciton
      (*distSgnd)[i][j] = signed_distance(x[i], y[j]);
    }
  }
}

template <size_t n, size_t m>
__always_inline double max_abs_diff(const unique_ptr<Array2D<n, m>> &mat1,
                                    const unique_ptr<Array2D<n, m>> &mat2) {
  auto res = 0.;

  for (size_t i = 0; i < x.size(); ++i) {
    for (size_t j = 0; j < y.size(); ++j) {
      res = max(res, abs((*mat1)[i + 3][j + 3] - (*mat2)[i + 3][j + 3]));
    }
  }

  return res;
}

int main(int argc, char *argv[]) {

  initialize();

  array<double, 4> mpd = {0, 0};

  auto w = make_unique<Array2D<x.size() + 6, y.size() + 6>>();

  for (size_t i = 0; i < x.size(); ++i) {
    for (size_t j = 0; j < y.size(); ++j) {
      // set initial approximation
      (*w)[i + 3][j + 3] =
          (__builtin_sin((x[i] + y[j] - .5) * M_PI * 2.)) * 1e-7 +
          (*distSgnd)[i][j];

      // Find maximum of derivarives
      mpd[0] = max(mpd[0], abs((*hx)[i][j] + (*mx)[i][j]));
      mpd[1] = max(mpd[1], abs((*hy)[i][j] + (*my)[i][j]));
    }
  }

  const auto dt = min(0.1 / max(mpd[0] / dx + mpd[1] / dy, L),
                      __builtin_pow(__builtin_sqrt(dx * dx + dy * dy), 1.5));
  // const auto dt =
  //     0.2 / max(mpd[0] / dx + mpd[1] / dy + mpd[2] / dy1 + mpd[3] / dy2,
  //     L);

  matrixContinuation(w);
  auto wnew0 = make_unique<Array2D<x.size() + 6, y.size() + 6>>();
  auto wnew1 = make_unique<Array2D<x.size() + 6, y.size() + 6>>();

  constexpr auto errorTolerance = 1e-8;
  auto currentError = errorTolerance;
  vector<double> err;
  auto ncounter = 0;
  auto start = chrono::high_resolution_clock::now();

  while (currentError >= errorTolerance) {

    // Heun's predictor-corrector method 2nd order
    performTimeStep(wnew0, w, dt);
    performTimeStep(wnew1, wnew0, dt);
    weightedSum(wnew0, w, wnew1, 1. / 2., 1. / 2.);

    // Heun's predictor-corrector method 3rd order
    // performTimeStep(wnew0, w, dt);
    // performTimeStep(wnew1, wnew0, dt);
    // weightedSum(wnew0, w, wnew1, 3. / 4., 1. / 4.);
    // performTimeStep(wnew1, wnew0, dt);
    // weightedSum(wnew0, w, wnew1, 1. / 3., 2. / 3.);

    currentError = max_abs_diff(w, wnew0);
    if (!err.empty() && currentError > err.back()) {
      cout << "AAAA" << "\n";
      cout << "Steps: " << ncounter << "\n";
      cout << "Error last: " << err.back() << " Current error: " << currentError
           << "\n";
    }
    err.push_back(currentError);
    w.swap(wnew0);
    ncounter++;
    if (ncounter % 5 == 0) {
      cout << "Steps: " << ncounter << "\n";
      cout << "Error: " << err.back() << "\n";
      cout << "Time: "
           << 1e-9 * chrono::duration_cast<chrono::nanoseconds>(
                         chrono::high_resolution_clock::now() - start)
                         .count()
           << "\n";
    }
  }
  auto end = chrono::high_resolution_clock::now();
  double time_taken =
      1e-9 * chrono::duration_cast<chrono::nanoseconds>(end - start).count();
  cout << endl;
  cout << "Final error: " << err.back() << endl;
  cout << "Steps: " << ncounter << endl;
  cout << "Main calculation took: " << fixed << time_taken << setprecision(9)
       << " sec" << endl;

  // const auto now =
  //     std::chrono::year_month_day(std::chrono::system_clock::now());
  // filesystem::path outputFile = "solution/finalSolution.csv";;
  filesystem::path outputFile =
      format("solution/finalSolution-{}.csv", end).c_str();
  filesystem::create_directories(outputFile.parent_path());
  ofstream finalSolution(outputFile);
  printSolutionExtended((*w), finalSolution);

  // system(format("cat {}", outputFile.c_str()).c_str());
  // ifstream finalSolutionReRead(outputFile);
  // cout << string((istreambuf_iterator<char>(finalSolutionReRead)),
  //                istreambuf_iterator<char>());
  return 0;
}
