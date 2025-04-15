#include <algorithm>
#include <array>
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

template <typename T, size_t n1, size_t n2, size_t n3, size_t n4>
using Array4DGeneric = array<array<array<array<T, n4>, n3>, n2>, n1>;

template <size_t n1, size_t n2, size_t n3, size_t n4>
// using Array4D = array<array<array<array<double, n4>, n3>, n2>, n1>;
using Array4D = Array4DGeneric<double, n1, n2, n3, n4>;

// parameters from Goldberg, Agusto 2021
// human
// constexpr auto N1 = 1e6;
// constexpr auto N2 = 1e7;
// constexpr auto Gamma1 = 1. / 14.;
// constexpr auto Gamma2 = 1. / 14.;
// constexpr auto BetaHV = 0.5;
// constexpr auto p11 = 0.5;
// constexpr auto p12 = 1 - p11;
// constexpr auto p22 = 0.5;
// constexpr auto p21 = 1 - p22;
constexpr auto N1 = 9.38e6;
constexpr auto N2 = 4.47e6;
constexpr auto Gamma1 = 0.063;
constexpr auto Gamma2 = 0.058;
constexpr auto p11 = 0.75;
constexpr auto p12 = 1 - p11;
constexpr auto p22 = 0.68;
constexpr auto p21 = 1 - p22;

constexpr auto H1 = p11 * N1 + p21 * N2;
constexpr auto H2 = p21 * N1 + p22 * N2;

// mosquito
// constexpr auto M1 = 1e7;
// constexpr auto M2 = 1e8;
// constexpr auto mu1 = 1. / 20.;
// constexpr auto mu2 = 1. / 20.;
// constexpr auto a1 = 1.;
// constexpr auto a2 = 1.;
constexpr auto M1 = 1.73e7;
constexpr auto M2 = 3e7;
constexpr auto mu1 = 0.032;
constexpr auto mu2 = 0.047;
constexpr auto a1 = 0.158;
constexpr auto a2 = 0.160;

constexpr auto BetaHV = 0.5;
constexpr auto BetaVH = 0.1;
constexpr auto tau = 10;

// multipliers in normalized model
const auto b11 = BetaVH * p11 * __builtin_exp(-mu1 * tau) * a1 * M1 / H1;
const auto b12 = BetaVH * p12 * __builtin_exp(-mu2 * tau) * a2 * M2 / H2;
const auto b21 = BetaVH * p21 * __builtin_exp(-mu1 * tau) * a1 * M1 / H1;
const auto b22 = BetaVH * p22 * __builtin_exp(-mu2 * tau) * a2 * M2 / H2;

constexpr auto c11 = BetaHV * p11 * a1 * N1 / H1;
constexpr auto c12 = BetaHV * p21 * a1 * N2 / H1;
constexpr auto c21 = BetaHV * p12 * a2 * N1 / H2;
constexpr auto c22 = BetaHV * p22 * a2 * N2 / H2;

// maximum infected
constexpr auto Imax1 = .15;
constexpr auto Imax2 = .15;
// set to reduce the size of the domain where the viability kernel lies
constexpr auto ymax1 = .1;
constexpr auto ymax2 = .1;

// control
// constexpr auto k = .6;
// constexpr auto umax1 = .3;
// constexpr auto umax2 = .4;
constexpr auto k = .4;
constexpr auto umax1 = .20;
constexpr auto umax2 = .14;

// threshold for I at maximum control
// const auto Is =
//     (alpha * Beta * nu * (1 - k * umax) * (1 - k * umax) - Gamma * mu) /
//     (alpha * Beta * nu * (1 - k * umax) * (1 - k * umax) +
//      Beta * Gamma * (1 - k * umax));

// lambda - lipschitz constant for the system
const auto Lexact = 3. * (b11 + b12 + b21 + b22 + c11 + c12 + c21 + c22) +
                    Gamma1 + Gamma2 + mu1 + mu2;
const auto L = 1.05 * Lexact;

// numerical space step
constexpr auto Nx1 = 25;
constexpr auto dx1 = Imax1 / (Nx1 - 1);
constexpr auto Nx2 = 25;
constexpr auto dx2 = Imax2 / (Nx2 - 1);
constexpr auto Ny1 = 25;
constexpr auto dy1 = ymax1 / (Ny1 - 1);
constexpr auto Ny2 = 25;
constexpr auto dy2 = ymax2 / (Ny2 - 1);

// x and y
constexpr auto padX1 = 4;
array<double, Nx1 + 2 * padX1> x1;
constexpr auto padX2 = 4;
array<double, Nx2 + 2 * padX2> x2;
constexpr auto padY1 = 4;
array<double, Ny1 + 2 * padY1> y1;
constexpr auto padY2 = 4;
array<double, Ny2 + 2 * padY2> y2;
unique_ptr<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>> hx1(nullptr);
unique_ptr<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>> hx2(nullptr);
unique_ptr<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>> hy1(nullptr);
unique_ptr<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>> hy2(nullptr);
unique_ptr<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>> mx1(nullptr);
unique_ptr<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>> mx2(nullptr);
unique_ptr<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>> my1(nullptr);
unique_ptr<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>> my2(nullptr);
unique_ptr<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>> wwx1(nullptr);
unique_ptr<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>> wwx2(nullptr);
unique_ptr<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>> wwy1(nullptr);
unique_ptr<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>> wwy2(nullptr);
unique_ptr<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>> ax1(nullptr);
unique_ptr<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>> ax2(nullptr);
unique_ptr<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>> ay1(nullptr);
unique_ptr<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>> ay2(nullptr);
unique_ptr<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>>
    distSgnd(nullptr);

unique_ptr<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>> Fx1p(nullptr);
unique_ptr<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>> Fx1m(nullptr);
unique_ptr<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>> Fx2p(nullptr);
unique_ptr<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>> Fx2m(nullptr);
unique_ptr<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>> Fy1p(nullptr);
unique_ptr<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>> Fy1m(nullptr);
unique_ptr<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>> Fy2p(nullptr);
unique_ptr<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>> Fy2m(nullptr);

unique_ptr<Array4D<x1.size() + 5, x2.size(), y1.size(), y2.size()>>
    wx1(nullptr);
unique_ptr<Array4D<x1.size(), x2.size() + 5, y1.size(), y2.size()>>
    wx2(nullptr);
unique_ptr<Array4D<x1.size(), x2.size(), y1.size() + 5, y2.size()>>
    wy1(nullptr);
unique_ptr<Array4D<x1.size(), x2.size(), y1.size(), y2.size() + 5>>
    wy2(nullptr);

unique_ptr<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>>
    hamiltonian(nullptr);

template <size_t n1, size_t n2, size_t n3, size_t n4>
__always_inline void
weightedSum(const unique_ptr<Array4D<n1, n2, n3, n4>> &res,
            const unique_ptr<Array4D<n1, n2, n3, n4>> &mat1,
            const unique_ptr<Array4D<n1, n2, n3, n4>> &mat2, const double &w1,
            const double &w2) noexcept {
  for (size_t i1 = 0; i1 < n1; ++i1) {
    for (size_t i2 = 0; i2 < n2; ++i2) {
      for (size_t i3 = 0; i3 < n3; ++i3) {
        for (size_t i4 = 0; i4 < n4; ++i4) {
          (*res)[i1][i2][i3][i4] =
              w1 * (*mat1)[i1][i2][i3][i4] + w2 * (*mat2)[i1][i2][i3][i4];
        }
      }
    }
  }
}

template <size_t n1, size_t n2, size_t n3, size_t n4>
__always_inline void printSolutionOmega(const Array4D<n1, n2, n3, n4> &arr,
                                        ostream &stream = std::cout,
                                        const string &delim = ", ") noexcept {
  for (size_t i1 = 0; i1 < Nx1; ++i1) {
    for (size_t i2 = 0; i2 < Nx2; ++i2) {
      for (size_t j1 = 0; j1 < Ny1; ++j1) {
        for (size_t j2 = 0; j2 < Ny2; ++j2) {
          stream << fixed << setprecision(8);
          stream << x1[i1 + padX1] << delim << x2[i2 + padX2] << delim
                 << y1[j1 + padY1] << delim << y2[j2 + padY2] << delim;
          stream << scientific;
          stream << arr[i1 + padX1 + 3][i2 + padX2 + 3][j1 + padY1 + 3]
                       [j2 + padY2 + 3]
                 << "\n";
        }
      }
    }
  }
}

template <size_t n1, size_t n2, size_t n3, size_t n4>
__always_inline void
printSolutionExtended(const Array4D<n1, n2, n3, n4> &arr,
                      ostream &stream = std::cout,
                      const string &delim = ", ") noexcept {
  for (size_t i1 = 0; i1 < x1.size(); ++i1) {
    for (size_t i2 = 0; i2 < x2.size(); ++i2) {
      for (size_t j1 = 0; j1 < y1.size(); ++j1) {
        for (size_t j2 = 0; j2 < y2.size(); ++j2) {
          stream << fixed << setprecision(8);
          stream << x1[i1] << delim << x2[i2] << delim << y1[j1] << delim
                 << y2[j2] << delim;
          stream << scientific;
          stream << arr[i1 + 3][i2 + 3][j1 + 3][j2 + 3] << "\n";
        }
      }
    }
  }
}

template <size_t n1, size_t n2, size_t n3, size_t n4>
__always_inline void
matrixContinuation(const unique_ptr<Array4D<n1, n2, n3, n4>> &mat) noexcept {

  for (size_t i1 = 0; i1 < n1 - 3; ++i1) {
    for (size_t i2 = 0; i2 < n2 - 3; ++i2) {
      for (size_t i3 = 0; i3 < n3 - 3; ++i3) {
        (*mat)[i1 + 3][i2 + 3][i3 + 3][0] = (*mat)[i1 + 3][i2 + 3][i3 + 3][1] =
            (*mat)[i1 + 3][i2 + 3][i3 + 3][2] =
                (*mat)[i1 + 3][i2 + 3][i3 + 3][3];
        (*mat)[i1 + 3][i2 + 3][i3 + 3][n4 - 1] =
            (*mat)[i1 + 3][i2 + 3][i3 + 3][n4 - 2] =
                (*mat)[i1 + 3][i2 + 3][i3 + 3][n4 - 3] =
                    (*mat)[i1 + 3][i2 + 3][i3 + 3][n4 - 4];
      }
    }
  }

  for (size_t i1 = 0; i1 < n1 - 3; ++i1) {
    for (size_t i2 = 0; i2 < n2 - 3; ++i2) {
      for (size_t i4 = 0; i4 < n4 - 3; ++i4) {
        (*mat)[i1 + 3][i2 + 3][0][i4 + 3] = (*mat)[i1 + 3][i2 + 3][1][i4 + 3] =
            (*mat)[i1 + 3][i2 + 3][2][i4 + 3] =
                (*mat)[i1 + 3][i2 + 3][3][i4 + 3];
        (*mat)[i1 + 3][i2 + 3][n3 - 1][i4 + 3] =
            (*mat)[i1 + 3][i2 + 3][n3 - 2][i4 + 3] =
                (*mat)[i1 + 3][i2 + 3][n3 - 3][i4 + 3] =
                    (*mat)[i1 + 3][i2 + 3][n3 - 4][i4 + 3];
      }
    }
  }

  for (size_t i1 = 0; i1 < n1 - 3; ++i1) {
    for (size_t i3 = 0; i3 < n3 - 3; ++i3) {
      for (size_t i4 = 0; i4 < n4 - 3; ++i4) {
        (*mat)[i1 + 3][0][i3 + 3][i4 + 3] = (*mat)[i1 + 3][1][i3 + 3][i4 + 3] =
            (*mat)[i1 + 3][2][i3 + 3][i4 + 3] =
                (*mat)[i1 + 3][3][i3 + 3][i4 + 3];
        (*mat)[i1 + 3][n2 - 1][i3 + 3][i4 + 3] =
            (*mat)[i1 + 3][n2 - 2][i3 + 3][i4 + 3] =
                (*mat)[i1 + 3][n2 - 3][i3 + 3][i4 + 3] =
                    (*mat)[i1 + 3][n2 - 4][i3 + 3][i4 + 3];
      }
    }
  }

  for (size_t i2 = 0; i2 < n2 - 3; ++i2) {
    for (size_t i3 = 0; i3 < n3 - 3; ++i3) {
      for (size_t i4 = 0; i4 < n4 - 3; ++i4) {
        (*mat)[0][i2 + 3][i3 + 3][i4 + 3] = (*mat)[1][i2 + 3][i3 + 3][i4 + 3] =
            (*mat)[2][i2 + 3][i3 + 3][i4 + 3] =
                (*mat)[3][i2 + 3][i3 + 3][i4 + 3];
        (*mat)[n1 - 1][i2 + 3][i3 + 3][i4 + 3] =
            (*mat)[n1 - 2][i2 + 3][i3 + 3][i4 + 3] =
                (*mat)[n1 - 3][i2 + 3][i3 + 3][i4 + 3] =
                    (*mat)[n1 - 4][i2 + 3][i3 + 3][i4 + 3];
      }
    }
  }

  for (size_t i1 = 0; i1 < 3; ++i1) {
    for (size_t i2 = 0; i2 < 3; ++i2) {
      for (size_t i3 = 0; i3 < 3; ++i3) {
        for (size_t i4 = 0; i4 < 3; ++i4) {
          (*mat)[i1][i2][i3][i4] = (*mat)[3][3][3][3];

          (*mat)[n1 - 3 + i1][i2][i3][i4] = (*mat)[n1 - 4][3][3][3];

          (*mat)[i1][n2 - 3 + i2][i3][i4] = (*mat)[3][n2 - 4][3][3];
          (*mat)[n1 - 3 + i1][n2 - 3 + i2][i3][i4] =
              (*mat)[n1 - 4][n2 - 4][3][3];

          (*mat)[i1][i2][n3 - 3 + i3][i4] = (*mat)[3][3][n3 - 4][3];
          (*mat)[n1 - 3 + i1][i2][n3 - 3 + i3][i4] =
              (*mat)[n1 - 4][3][n3 - 4][3];
          (*mat)[i1][n2 - 3 + i2][n3 - 3 + i3][i4] =
              (*mat)[3][n2 - 4][n3 - 4][3];
          (*mat)[n1 - 3 + i1][n2 - 3 + i2][n3 - 3 + i3][i4] =
              (*mat)[n1 - 4][n2 - 4][n3 - 4][3];

          (*mat)[i1][i2][i3][n4 - 3 + i4] = (*mat)[3][3][3][n4 - 4];
          (*mat)[n1 - 3 + i1][i2][i3][n4 - 3 + i4] =
              (*mat)[n1 - 4][3][3][n4 - 4];
          (*mat)[i1][n2 - 3 + i2][i3][n4 - 3 + i4] =
              (*mat)[3][n2 - 4][3][n4 - 4];
          (*mat)[n1 - 3 + i1][n2 - 3 + i2][i3][n4 - 3 + i4] =
              (*mat)[n1 - 4][n2 - 4][3][n4 - 4];
          (*mat)[i1][i2][n3 - 3 + i3][n4 - 3 + i4] =
              (*mat)[3][3][n3 - 4][n4 - 4];
          (*mat)[n1 - 3 + i1][i2][n3 - 3 + i3][n4 - 3 + i4] =
              (*mat)[n1 - 4][3][n3 - 4][n4 - 4];
          (*mat)[i1][n2 - 3 + i2][n3 - 3 + i3][n4 - 3 + i4] =
              (*mat)[3][n2 - 4][n3 - 4][n4 - 4];
          (*mat)[n1 - 3 + i1][n2 - 3 + i2][n3 - 3 + i3][n4 - 3 + i4] =
              (*mat)[n1 - 4][n2 - 4][n3 - 4][n4 - 4];
        }
      }
    }
  }
}

__always_inline double weno(const double &v1, const double &v2,
                            const double &v3, const double &v4,
                            const double &v5) noexcept {
  // smoothness of stencils
  const auto S1 = 13. / 12. * (v1 - 2. * v2 + v3) * (v1 - 2. * v2 + v3) +
                  1. / 4. * (v1 - 4. * v2 + 3. * v3) * (v1 - 4. * v2 + 3. * v3);
  const auto S2 = 13. / 12. * (v2 - 2. * v3 + v4) * (v2 - 2. * v3 + v4) +
                  1. / 4. * (v2 - v4) * (v2 - v4);
  const auto S3 = 13. / 12. * (v3 - 2. * v4 + v5) * (v3 - 2. * v4 + v5) +
                  1. / 4. * (3. * v3 - 4. * v4 + v5) * (3. * v3 - 4. * v4 + v5);
  // find the maximum of vi's
  const auto maxv2 = max({v1 * v1, v2 * v2, v3 * v3, v4 * v4, v5 * v5});
  const auto eps1 = 1e-6 * maxv2 + 1e-99;
  // weights
  const auto a1 = .1 / ((S1 + eps1) * (S1 + eps1));
  const auto a2 = .6 / ((S2 + eps1) * (S2 + eps1));
  const auto a3 = .3 / ((S3 + eps1) * (S3 + eps1));
  const auto w1 = a1 / (a1 + a2 + a3);
  const auto w2 = a2 / (a1 + a2 + a3);
  const auto w3 = a3 / (a1 + a2 + a3);
  const auto f1x = v1 / 3. - 7. * v2 / 6 + 11. * v3 / 6.;
  const auto f2x = -v2 / 6. + 5. * v3 / 6. + v4 / 3.;
  const auto f3x = v3 / 3. + 5. * v4 / 6. - v5 / 6.;
  return w1 * f1x + w2 * f2x + w3 * f3x;
}

template <size_t n1, size_t n2, size_t n3, size_t n4>
__always_inline void approximate_derivatives(
    const unique_ptr<Array4D<n1, n2, n3, n4>> &mat) noexcept {
  for (size_t i1 = 0; i1 < n1 - 1; ++i1) {
    for (size_t i2 = 0; i2 < n2 - 6; ++i2) {
      for (size_t i3 = 0; i3 < n3 - 6; ++i3) {
        for (size_t i4 = 0; i4 < n4 - 6; ++i4) {
          (*wx1)[i1][i2][i3][i4] = ((*mat)[i1 + 1][i2 + 3][i3 + 3][i4 + 3] -
                                    (*mat)[i1][i2 + 3][i3 + 3][i4 + 3]) /
                                   dx1;
        }
      }
    }
  }
  for (size_t i1 = 0; i1 < n1 - 6; ++i1) {
    for (size_t i2 = 0; i2 < n2 - 1; ++i2) {
      for (size_t i3 = 0; i3 < n3 - 6; ++i3) {
        for (size_t i4 = 0; i4 < n4 - 6; ++i4) {
          (*wx2)[i1][i2][i3][i4] = ((*mat)[i1 + 3][i2 + 1][i3 + 3][i4 + 3] -
                                    (*mat)[i1 + 3][i2][i3 + 3][i4 + 3]) /
                                   dx2;
        }
      }
    }
  }
  for (size_t i1 = 0; i1 < n1 - 6; ++i1) {
    for (size_t i2 = 0; i2 < n2 - 6; ++i2) {
      for (size_t i3 = 0; i3 < n3 - 1; ++i3) {
        for (size_t i4 = 0; i4 < n4 - 6; ++i4) {
          (*wy1)[i1][i2][i3][i4] = ((*mat)[i1 + 3][i2 + 3][i3 + 1][i4 + 3] -
                                    (*mat)[i1 + 3][i2 + 3][i3][i4 + 3]) /
                                   dy1;
        }
      }
    }
  }
  for (size_t i1 = 0; i1 < n1 - 6; ++i1) {
    for (size_t i2 = 0; i2 < n2 - 6; ++i2) {
      for (size_t i3 = 0; i3 < n3 - 6; ++i3) {
        for (size_t i4 = 0; i4 < n4 - 1; ++i4) {
          (*wy2)[i1][i2][i3][i4] = ((*mat)[i1 + 3][i2 + 3][i3 + 3][i4 + 1] -
                                    (*mat)[i1 + 3][i2 + 3][i3 + 3][i4]) /
                                   dy2;
        }
      }
    }
  }

  for (size_t i1 = 0; i1 < n1 - 6; ++i1) {
    for (size_t i2 = 0; i2 < n2 - 6; ++i2) {
      for (size_t i3 = 0; i3 < n3 - 6; ++i3) {
        for (size_t i4 = 0; i4 < n4 - 6; ++i4) {
          (*Fx1m)[i1][i2][i3][i4] =
              weno((*wx1)[i1][i2][i3][i4], (*wx1)[i1 + 1][i2][i3][i4],
                   (*wx1)[i1 + 2][i2][i3][i4], (*wx1)[i1 + 3][i2][i3][i4],
                   (*wx1)[i1 + 4][i2][i3][i4]);
          (*Fx1p)[i1][i2][i3][i4] =
              weno((*wx1)[i1 + 5][i2][i3][i4], (*wx1)[i1 + 4][i2][i3][i4],
                   (*wx1)[i1 + 3][i2][i3][i4], (*wx1)[i1 + 2][i2][i3][i4],
                   (*wx1)[i1 + 1][i2][i3][i4]);

          (*Fx2m)[i1][i2][i3][i4] =
              weno((*wx2)[i1][i2][i3][i4], (*wx2)[i1][i2 + 1][i3][i4],
                   (*wx2)[i1][i2 + 2][i3][i4], (*wx2)[i1][i2 + 3][i3][i4],
                   (*wx2)[i1][i2 + 4][i3][i4]);
          (*Fx2p)[i1][i2][i3][i4] =
              weno((*wx2)[i1][i2 + 5][i3][i4], (*wx2)[i1][i2 + 4][i3][i4],
                   (*wx2)[i1][i2 + 3][i3][i4], (*wx2)[i1][i2 + 2][i3][i4],
                   (*wx2)[i1][i2 + 1][i3][i4]);

          (*Fy1m)[i1][i2][i3][i4] =
              weno((*wy1)[i1][i2][i3][i4], (*wy1)[i1][i2][i3 + 1][i4],
                   (*wy1)[i1][i2][i3 + 2][i4], (*wy1)[i1][i2][i3 + 3][i4],
                   (*wy1)[i1][i2][i3 + 4][i4]);
          (*Fy1p)[i1][i2][i3][i4] =
              weno((*wy1)[i1][i2][i3 + 5][i4], (*wy1)[i1][i2][i3 + 4][i4],
                   (*wy1)[i1][i2][i3 + 3][i4], (*wy1)[i1][i2][i3 + 2][i4],
                   (*wy1)[i1][i2][i3 + 1][i4]);

          (*Fy2m)[i1][i2][i3][i4] =
              weno((*wy2)[i1][i2][i3][i4], (*wy2)[i1][i2][i3][i4 + 1],
                   (*wy2)[i1][i2][i3][i4 + 2], (*wy2)[i1][i2][i3][i4 + 3],
                   (*wy2)[i1][i2][i3][i4 + 4]);
          (*Fy2p)[i1][i2][i3][i4] =
              weno((*wy2)[i1][i2][i3][i4 + 5], (*wy2)[i1][i2][i3][i4 + 4],
                   (*wy2)[i1][i2][i3][i4 + 3], (*wy2)[i1][i2][i3][i4 + 2],
                   (*wy2)[i1][i2][i3][i4 + 1]);
        }
      }
    }
  }
}

__always_inline double hamiltonianInDirection(const double &fp,
                                              const double &fm,
                                              const double &ww, const double &h,
                                              const double &m,
                                              const double &a) noexcept {
  double res = 0;
  const auto changeInSign = ww * (fp * fm < 1e-99);
  const auto interm = h + (m * (fp + fm) / 2. > 1e-99) * m;
  res += (1 - changeInSign) *
         ((interm < 1e-99) * interm * fp + (interm > 1e-99) * interm * fm);
  res += changeInSign * ((interm - a) * fp + (interm + a) * fm) / 2.;
  return res;
}

__always_inline void calcHamiltonian() noexcept {

  for (size_t i1 = 0; i1 < x1.size(); ++i1) {
    for (size_t i2 = 0; i2 < x2.size(); ++i2) {
      for (size_t i3 = 0; i3 < y1.size(); ++i3) {
        for (size_t i4 = 0; i4 < y2.size(); ++i4) {
          auto hamilt = 0.;

          hamilt += hamiltonianInDirection(
              (*Fx1p)[i1][i2][i3][i4], (*Fx1m)[i1][i2][i3][i4],
              (*wwx1)[i1][i2][i3][i4], (*hx1)[i1][i2][i3][i4],
              (*mx1)[i1][i2][i3][i4], (*ax1)[i1][i2][i3][i4]);
          hamilt += hamiltonianInDirection(
              (*Fx2p)[i1][i2][i3][i4], (*Fx2m)[i1][i2][i3][i4],
              (*wwx2)[i1][i2][i3][i4], (*hx2)[i1][i2][i3][i4],
              (*mx2)[i1][i2][i3][i4], (*ax2)[i1][i2][i3][i4]);
          hamilt += hamiltonianInDirection(
              (*Fy1p)[i1][i2][i3][i4], (*Fy1m)[i1][i2][i3][i4],
              (*wwy1)[i1][i2][i3][i4], (*hy1)[i1][i2][i3][i4],
              (*my1)[i1][i2][i3][i4], (*ay1)[i1][i2][i3][i4]);
          hamilt += hamiltonianInDirection(
              (*Fy2p)[i1][i2][i3][i4], (*Fy2m)[i1][i2][i3][i4],
              (*wwy2)[i1][i2][i3][i4], (*hy2)[i1][i2][i3][i4],
              (*my2)[i1][i2][i3][i4], (*ay2)[i1][i2][i3][i4]);

          (*hamiltonian)[i1][i2][i3][i4] = hamilt;
        }
      }
    }
  }
}

template <size_t n1, size_t n2, size_t n3, size_t n4>
__always_inline void
performTimeStep(const unique_ptr<Array4D<n1, n2, n3, n4>> &tp1Matrix,
                const unique_ptr<Array4D<n1, n2, n3, n4>> &tMatrix,
                const double &dt) noexcept {
  approximate_derivatives(tMatrix);
  calcHamiltonian();

  for (size_t i1 = 0; i1 < n1 - 6; ++i1) {
    for (size_t i2 = 0; i2 < n2 - 6; ++i2) {
      for (size_t j1 = 0; j1 < n3 - 6; ++j1) {
        for (size_t j2 = 0; j2 < n4 - 6; ++j2) {
          (*tp1Matrix)[i1 + 3][i2 + 3][j1 + 3][j2 + 3] =
              max((1 - L * dt) * (*tMatrix)[i1 + 3][i2 + 3][j1 + 3][j2 + 3] -
                      dt * (*hamiltonian)[i1][i2][j1][j2],
                  (*distSgnd)[i1][i2][j1][j2]);
        }
      }
    }
  }

  matrixContinuation(tp1Matrix);
}

__always_inline double signed_distance(const double &x1, const double &x2) {
  if (x1 >= 0 && x2 >= 0) {
    if (x1 <= Imax1) {
      if (x2 <= Imax2) {
        return -min({abs(x1 - Imax1), abs(x2 - Imax2)});
        // return -min({abs(x1 - Imax1), abs(x2 - Imax2), abs(x1), abs(x2)});
      } else {
        return abs(x2 - Imax2);
      }
    } else if (x2 <= Imax2) {
      return abs(x1 - Imax1);
    } else
      return __builtin_sqrt((x1 - Imax1) * (x1 - Imax1) +
                            (x2 - Imax2) * (x2 - Imax2));
  } else if (x1 < 0) {
    if (x2 >= 0) {
      if (x2 <= Imax2) {
        return abs(x1);
      } else
        return __builtin_sqrt(x1 * x1 + (x2 - Imax2) * (x2 - Imax2));
    } else
      return __builtin_sqrt(x1 * x1 + x2 * x2);
  } else {
    if (x1 <= Imax1) {
      return abs(x2);
    } else
      return __builtin_sqrt((x1 - Imax1) * (x1 - Imax1) + x2 * x2);
  }
}

__always_inline void initialize() {
  hx1 = make_unique<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>>();
  hx2 = make_unique<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>>();
  hy1 = make_unique<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>>();
  hy2 = make_unique<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>>();
  mx1 = make_unique<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>>();
  mx2 = make_unique<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>>();
  my1 = make_unique<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>>();
  my2 = make_unique<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>>();
  wwx1 = make_unique<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>>();
  wwx2 = make_unique<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>>();
  wwy1 = make_unique<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>>();
  wwy2 = make_unique<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>>();
  ax1 = make_unique<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>>();
  ax2 = make_unique<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>>();
  ay1 = make_unique<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>>();
  ay2 = make_unique<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>>();

  Fx1p = make_unique<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>>();
  Fx1m = make_unique<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>>();
  Fx2p = make_unique<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>>();
  Fx2m = make_unique<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>>();
  Fy1p = make_unique<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>>();
  Fy1m = make_unique<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>>();
  Fy2p = make_unique<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>>();
  Fy2m = make_unique<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>>();

  hamiltonian =
      make_unique<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>>();

  distSgnd = make_unique<Array4D<x1.size(), x2.size(), y1.size(), y2.size()>>();

  wx1 = make_unique<Array4D<x1.size() + 5, x2.size(), y1.size(), y2.size()>>();
  wx2 = make_unique<Array4D<x1.size(), x2.size() + 5, y1.size(), y2.size()>>();
  wy1 = make_unique<Array4D<x1.size(), x2.size(), y1.size() + 5, y2.size()>>();
  wy2 = make_unique<Array4D<x1.size(), x2.size(), y1.size(), y2.size() + 5>>();

  generate(x1.begin(), x1.end(),
           [n = 0]() mutable { return (-padX1 + n++) * dx1; });
  x1[padX1] = 0;
  x1[padX1 + Nx1 - 1] = Imax1;
  generate(x2.begin(), x2.end(),
           [n = 0]() mutable { return (-padX2 + n++) * dx2; });
  x2[padX2] = 0;
  x2[padX2 + Nx2 - 1] = Imax2;
  generate(y1.begin(), y1.end(),
           [n = 0]() mutable { return (-padY1 + n++) * dy1; });
  y1[padY1] = 0;
  y1[padY1 + Ny1 - 1] = ymax1;
  generate(y2.begin(), y2.end(),
           [n = 0]() mutable { return (-padY2 + n++) * dy2; });
  y2[padY2] = 0;
  y2[padY2 + Ny2 - 1] = ymax2;

  // generate(x1.begin(), x1.end(), [n = 0]() mutable { return (n++) * dx1;
  // }); x1[Nx1 - 1] = Imax1; generate(x2.begin(), x2.end(), [n = 0]() mutable
  // { return (n++) * dx2; }); x2[Nx2 - 1] = Imax2; generate(y1.begin(),
  // y1.end(), [n = 0]() mutable { return (n++) * dy1; }); y1[Ny1 - 1] =
  // ymax1; generate(y2.begin(), y2.end(), [n = 0]() mutable { return (n++) *
  // dy2; }); y2[Ny2 - 1] = ymax2;
}

template <size_t n1, size_t n2, size_t n3, size_t n4>
__always_inline double
max_abs_diff(const unique_ptr<Array4D<n1, n2, n3, n4>> &mat1,
             const unique_ptr<Array4D<n1, n2, n3, n4>> &mat2) {
  auto res = 0.;
  // auto solSize = 0.;
  for (size_t i1 = 0; i1 < x1.size(); ++i1) {
    for (size_t i2 = 0; i2 < x2.size(); ++i2) {
      for (size_t j1 = 0; j1 < y1.size(); ++j1) {
        for (size_t j2 = 0; j2 < y2.size(); ++j2) {
          res = max(res, abs((*mat1)[i1 + 3][i2 + 3][j1 + 3][j2 + 3] -
                             (*mat2)[i1 + 3][i2 + 3][j1 + 3][j2 + 3]));
          // solSize += abs((*wnew0)[i1 + 3][i2 + 3][j1 + 3][j2 + 3]);
          // currentError = max(currentError,
          //                    abs((*w)[i1 + 3][i2 + 3][j1 + 3][j2 + 3] -
          //                        (*wnew0)[i1 + 3][i2 + 3][j1 + 3][j2 +
          //                        3]));

          // currentError = max(currentError,
          //                    abs(((*w)[i1 + 3][i2 + 3][j1 + 3][j2 + 3] -
          //                         (*wnew0)[i1 + 3][i2 + 3][j1 + 3][j2 + 3])
          //                         /
          //                        (*w)[i1 + 3][i2 + 3][j1 + 3][j2 + 3]));
        }
      }
    }
  }
  return res;
}

int main(int argc, char *argv[]) {

  initialize();

  array<double, 4> mpd = {0, 0, 0, 0};

  auto w = make_unique<
      Array4D<x1.size() + 6, x2.size() + 6, y1.size() + 6, y2.size() + 6>>();

  // compute the partial derivatives
  for (size_t i1 = 0; i1 < x1.size(); ++i1) {
    for (size_t i2 = 0; i2 < x2.size(); ++i2) {
      for (size_t j1 = 0; j1 < y1.size(); ++j1) {
        for (size_t j2 = 0; j2 < y2.size(); ++j2) {
          (*hx1)[i1][i2][j1][j2] =
              Gamma1 * x1[i1] - (1 - x1[i1]) * (b11 * y1[j1] + b12 * y2[j2]);
          (*hx2)[i1][i2][j1][j2] =
              Gamma2 * x2[i2] - (1 - x2[i2]) * (b21 * y1[j1] + b22 * y2[j2]);
          (*hy1)[i1][i2][j1][j2] =
              mu1 * y1[j1] - (1 - y1[j1]) * (c11 * x1[i1] + c12 * x2[i2]);
          (*hy2)[i1][i2][j1][j2] =
              mu2 * y2[j2] - (1 - y2[j2]) * (c21 * x1[i1] + c22 * x2[i2]);

          (*mx1)[i1][i2][j1][j2] =
              k * umax1 * (1 - x1[i1]) * (b11 * y1[j1] + b12 * y2[j2]);
          (*mx2)[i1][i2][j1][j2] =
              k * umax2 * (1 - x2[i2]) * (b21 * y1[j1] + b22 * y2[j2]);
          (*my1)[i1][i2][j1][j2] =
              k * (1 - y1[j1]) * (c11 * umax1 * x1[i1] + c12 * umax2 * x2[i2]);
          (*my2)[i1][i2][j1][j2] =
              k * (1 - y2[j2]) * (c21 * umax1 * x1[i1] + c22 * umax2 * x2[i2]);

          // find maximum range of gridpoints
          (*wwx1)[i1][i2][j1][j2] =
              ((*hx1)[i1][i2][j1][j2] *
               ((*hx1)[i1][i2][j1][j2] + (*mx1)[i1][i2][j1][j2])) < 1e-99;
          (*ax1)[i1][i2][j1][j2] =
              max(abs((*hx1)[i1][i2][j1][j2]),
                  abs((*hx1)[i1][i2][j1][j2] + (*mx1)[i1][i2][j1][j2]));
          (*wwx2)[i1][i2][j1][j2] =
              ((*hx2)[i1][i2][j1][j2] *
               ((*hx2)[i1][i2][j1][j2] + (*mx2)[i1][i2][j1][j2])) < 1e-99;
          (*ax2)[i1][i2][j1][j2] =
              max(abs((*hx2)[i1][i2][j1][j2]),
                  abs((*hx2)[i1][i2][j1][j2] + (*mx2)[i1][i2][j1][j2]));
          (*wwy1)[i1][i2][j1][j2] =
              ((*hy1)[i1][i2][j1][j2] *
               ((*hy1)[i1][i2][j1][j2] + (*my1)[i1][i2][j1][j2])) < 1e-99;
          (*ay1)[i1][i2][j1][j2] =
              max(abs((*hy1)[i1][i2][j1][j2]),
                  abs((*hy1)[i1][i2][j1][j2] + (*my1)[i1][i2][j1][j2]));
          (*wwy2)[i1][i2][j1][j2] =
              ((*hy2)[i1][i2][j1][j2] *
               ((*hy2)[i1][i2][j1][j2] + (*my2)[i1][i2][j1][j2])) < 1e-99;
          (*ay2)[i1][i2][j1][j2] =
              max(abs((*hy2)[i1][i2][j1][j2]),
                  abs((*hy2)[i1][i2][j1][j2] + (*my2)[i1][i2][j1][j2]));

          // set signed distance
          (*distSgnd)[i1][i2][j1][j2] = signed_distance(x1[i1], x2[i2]);

          // set initial approximation
          (*w)[i1 + 3][i2 + 3][j1 + 3][j2 + 3] =
              (__builtin_sin((x1[i1] + x2[i2] + y1[j1] + y2[j2] - .5) * M_PI *
                             2.)) *
                  1e-5 +
              (*distSgnd)[i1][i2][j1][j2];

          // Find maximum of derivarives
          mpd[0] =
              max(mpd[0], abs((*hx1)[i1][i2][j1][j2] + (*mx1)[i1][i2][j1][j2]));
          mpd[1] =
              max(mpd[1], abs((*hx2)[i1][i2][j1][j2] + (*mx2)[i1][i2][j1][j2]));
          mpd[2] =
              max(mpd[2], abs((*hy1)[i1][i2][j1][j2] + (*my1)[i1][i2][j1][j2]));
          mpd[3] =
              max(mpd[3], abs((*hy2)[i1][i2][j1][j2] + (*my2)[i1][i2][j1][j2]));
        }
      }
    }
  }

  const auto dt =
      0.9 / max(mpd[0] / dx1 + mpd[1] / dx2 + mpd[2] / dy1 + mpd[3] / dy2, L);

  matrixContinuation(w);
  auto wnew0 = make_unique<
      Array4D<x1.size() + 6, x2.size() + 6, y1.size() + 6, y2.size() + 6>>();
  auto wnew1 = make_unique<
      Array4D<x1.size() + 6, x2.size() + 6, y1.size() + 6, y2.size() + 6>>();

  constexpr auto errtol = 1e-10;
  auto currentError = errtol;
  vector<double> err;
  auto ncounter = 0;
  auto start = chrono::high_resolution_clock::now();

  while (currentError >= errtol) {

    // Heun's predictor-corrector method 2nd order
    // performTimeStep(wnew0, w, dt);
    // performTimeStep(wnew1, wnew0, dt);
    // weightedSum(wnew0, w, wnew1, 1. / 2., 1. / 2.);

    // Heun's predictor-corrector method 3rd order
    performTimeStep(wnew0, w, dt);
    performTimeStep(wnew1, wnew0, dt);
    weightedSum(wnew0, w, wnew1, 3. / 4., 1. / 4.);
    performTimeStep(wnew1, wnew0, dt);
    weightedSum(wnew0, w, wnew1, 1. / 3., 2. / 3.);

    currentError = max_abs_diff(w, wnew0);
    if (!err.empty() && currentError > err.back()) {
      cout << "AAAA" << "\n";
      cout << "Steps: " << ncounter << "\n";
      cout << "Error last: " << err.back() << "Current error: " << currentError
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
}
