#include <algorithm>
#include <array>
#include <chrono>
// #include <cmath>
#include <filesystem>
#include <format>
#include <fstream>
#include <iomanip>
#include <iostream>
// #include <iterator>
#include <tuple>
#include <vector>
using namespace std;

// #undef y1
constexpr double M_PI = 3.141592653589793115997963468544185161590576171875;

template <typename T, size_t n1, size_t n2, size_t n3, size_t n4>
using Array4DGeneric = array<array<array<array<T, n4>, n3>, n2>, n1>;

template <size_t n1, size_t n2, size_t n3, size_t n4>
// using Array4D = array<array<array<array<double, n4>, n3>, n2>, n1>;
using Array4D = Array4DGeneric<double, n1, n2, n3, n4>;

template <size_t n1, size_t n2, size_t n3, size_t n4>
struct derivativeApproximations {
  Array4D<n1, n2, n3, n4> fx1p;
  Array4D<n1, n2, n3, n4> fx1m;
  Array4D<n1, n2, n3, n4> fx2p;
  Array4D<n1, n2, n3, n4> fx2m;
  Array4D<n1, n2, n3, n4> fy1p;
  Array4D<n1, n2, n3, n4> fy1m;
  Array4D<n1, n2, n3, n4> fy2p;
  Array4D<n1, n2, n3, n4> fy2m;
};

// parameters from Goldberg, Agusto 2021
// human
constexpr auto N1 = 1e6;
constexpr auto N2 = 1e7;
constexpr auto Gamma1 = 1. / 14.;
constexpr auto Gamma2 = 1. / 14.;
constexpr auto BetaHV = 0.5;
constexpr auto p11 = 0.5;
constexpr auto p12 = 1 - p11;
constexpr auto p22 = 0.5;
constexpr auto p21 = 1 - p22;
constexpr auto H1 = p11 * N1 + p21 * N2;
constexpr auto H2 = p12 * N1 + p22 * N2;

// mosquito
constexpr auto M1 = 1e7;
constexpr auto M2 = 1e8;
constexpr auto mu1 = 1. / 20.;
constexpr auto mu2 = 1. / 20.;
constexpr auto a1 = 1.;
constexpr auto a2 = 1.;
constexpr auto BetaVH = 0.1;
constexpr auto tau = 10.;

const auto b11 = BetaVH * p11 * __builtin_exp(-mu1 * tau) * a1 * M1 / H1;
const auto b12 = BetaVH * p12 * __builtin_exp(-mu2 * tau) * a2 * M2 / H2;
const auto b21 = BetaVH * p21 * __builtin_exp(-mu1 * tau) * a1 * M1 / H1;
const auto b22 = BetaVH * p22 * __builtin_exp(-mu2 * tau) * a2 * M2 / H2;

constexpr auto c11 = BetaHV * p11 * a1 * N1 / H1;
constexpr auto c12 = BetaHV * p21 * a1 * N2 / H1;
constexpr auto c21 = BetaHV * p12 * a2 * N1 / H2;
constexpr auto c22 = BetaHV * p22 * a2 * N2 / H2;

// maximum infected
constexpr auto Imax1 = .1;
constexpr auto Imax2 = .1;
// set to reduce the size of the domain where the viability kernel lies
constexpr auto ymax1 = .5;
constexpr auto ymax2 = .5;

// control
constexpr auto k = .6;
constexpr auto umax1 = .3;
constexpr auto umax2 = .4;

// predefine
constexpr auto errtol = 1e-10;

// lambda - lipschitz constant for the system
const auto Lexact = 3. * (b11 + b12 + b21 + b22 + c11 + c12 + c21 + c22) +
                    Gamma1 + Gamma2 + mu1 + mu2;
const auto L = 1.05 * Lexact;

// maximum of partial derivatives, Hamiltonian
const array<double, 4> mpd = {
    1 + Gamma1 + (b11 + b12) * (1 + k * umax1),
    1 + Gamma2 + (b21 + b22) * (1 + k * umax2),
    1 + mu1 + (c11 * (1 + k * umax1) + c12 * (1 + k * umax2)),
    1 + mu2 + (c21 * (1 + k * umax1) + c22 * (1 + k * umax2))};
// mpd = mpd + 1; ????

// numerical space step
constexpr auto Nx1 = 12;
constexpr auto dx1 = Imax1 / (Nx1 - 1);
constexpr auto Nx2 = 12;
constexpr auto dx2 = Imax2 / (Nx2 - 1);
constexpr auto Ny1 = 8;
constexpr auto dy1 = ymax1 / (Ny1 - 1);
constexpr auto Ny2 = 8;
constexpr auto dy2 = ymax2 / (Ny2 - 1);

// numerical time step dt < 1/ll (!)
// numerical time step
// const auto dt = 0.9 * min({1. / (mpd[0] / dx + mpd[1] / dy), 1. / ll});
const auto dt =
    0.9 / max(mpd[0] / dx1 + mpd[1] / dx2 + mpd[2] / dy1 + mpd[3] / dy2, L);

// x and y
constexpr auto padX1 = 3;
array<double, Nx1> x1;
constexpr auto padX2 = 3;
array<double, Nx2> x2;
constexpr auto padY1 = 3;
array<double, Ny1> y1;
constexpr auto padY2 = 3;
array<double, Ny2> y2;
static Array4D<x1.size(), x2.size(), y1.size(), y2.size()> hx1;
static Array4D<x1.size(), x2.size(), y1.size(), y2.size()> hx2;
static Array4D<x1.size(), x2.size(), y1.size(), y2.size()> hy1;
static Array4D<x1.size(), x2.size(), y1.size(), y2.size()> hy2;
static Array4D<x1.size(), x2.size(), y1.size(), y2.size()> mx1;
static Array4D<x1.size(), x2.size(), y1.size(), y2.size()> mx2;
static Array4D<x1.size(), x2.size(), y1.size(), y2.size()> my1;
static Array4D<x1.size(), x2.size(), y1.size(), y2.size()> my2;
static Array4D<x1.size(), x2.size(), y1.size(), y2.size()> wwx1;
static Array4D<x1.size(), x2.size(), y1.size(), y2.size()> wwx2;
static Array4D<x1.size(), x2.size(), y1.size(), y2.size()> wwy1;
static Array4D<x1.size(), x2.size(), y1.size(), y2.size()> wwy2;
static Array4D<x1.size(), x2.size(), y1.size(), y2.size()> ax1;
static Array4D<x1.size(), x2.size(), y1.size(), y2.size()> ax2;
static Array4D<x1.size(), x2.size(), y1.size(), y2.size()> ay1;
static Array4D<x1.size(), x2.size(), y1.size(), y2.size()> ay2;
static Array4D<x1.size(), x2.size(), y1.size(), y2.size()> distSgnd;
Array4D<x1.size() + 2 * padX1, x2.size() + 2 * padX2, y1.size() + 2 * padY1,
        y2.size() + 2 * padY2>
    w;

template <size_t n1, size_t n2, size_t n3, size_t n4>
Array4D<n1, n2, n3, n4> halfSum(const Array4D<n1, n2, n3, n4> &mat1,
                                const Array4D<n1, n2, n3, n4> &mat2) noexcept {
  Array4D<n1, n2, n3, n4> res;
  for (size_t i1 = 0; i1 < n1; ++i1) {
    for (size_t i2 = 0; i2 < n2; ++i2) {
      for (size_t i3 = 0; i3 < n3; ++i3) {
        for (size_t i4 = 0; i4 < n4; ++i4) {
          res[i1][i2][i3][i4] =
              (mat1[i1][i2][i3][i4] + mat2[i1][i2][i3][i4]) / 2.;
        }
      }
    }
  }
  return res;
}

template <size_t n1, size_t n2, size_t n3, size_t n4>
void printSolution(const Array4D<n1, n2, n3, n4> &arr,
                   ostream &stream = std::cout,
                   const string &delim = ", ") noexcept {
  for (size_t i1 = 0; i1 < Nx1; ++i1) {
    for (size_t i2 = 0; i2 < Nx2; ++i2) {
      for (size_t j1 = 0; j1 < Ny1; ++j1) {
        for (size_t j2 = 0; j2 < Ny2; ++j2) {
          stream << x1[i1] << delim << x2[i2] << delim << y1[j1] << delim
                 << y2[j2] << delim
                 << arr[i1 + padX1][i2 + padX2][j1 + padY1][j2 + padY2] << "\n";
        }
      }
    }
  }
}

template <size_t n1, size_t n2, size_t n3, size_t n4>
void matrixContinuation(Array4D<n1, n2, n3, n4> &mat) noexcept {

  for (size_t i1 = 0; i1 < n1 - padX1; ++i1) {
    for (size_t i2 = 0; i2 < n2 - padX2; ++i2) {
      for (size_t i3 = 0; i3 < n3 - padY1; ++i3) {
        mat[i1 + 3][i2 + 3][i3 + 3][0] = mat[i1 + 3][i2 + 3][i3 + 3][1] =
            mat[i1 + 3][i2 + 3][i3 + 3][2] = mat[i1 + 3][i2 + 3][i3 + 3][3];
        mat[i1 + 3][i2 + 3][i3 + 3][n4 - 1] =
            mat[i1 + 3][i2 + 3][i3 + 3][n4 - 2] =
                mat[i1 + 3][i2 + 3][i3 + 3][n4 - 3] =
                    mat[i1 + 3][i2 + 3][i3 + 3][n4 - 4];
      }
    }
  }

  for (size_t i1 = 0; i1 < n1 - padX1; ++i1) {
    for (size_t i2 = 0; i2 < n2 - padX2; ++i2) {
      for (size_t i4 = 0; i4 < n4 - padY2; ++i4) {
        mat[i1 + 3][i2 + 3][0][i4 + 3] = mat[i1 + 3][i2 + 3][1][i4 + 3] =
            mat[i1 + 3][i2 + 3][2][i4 + 3] = mat[i1 + 3][i2 + 3][3][i4 + 3];
        mat[i1 + 3][i2 + 3][n3 - 1][i4 + 3] =
            mat[i1 + 3][i2 + 3][n3 - 2][i4 + 3] =
                mat[i1 + 3][i2 + 3][n3 - 3][i4 + 3] =
                    mat[i1 + 3][i2 + 3][n3 - 4][i4 + 3];
      }
    }
  }

  for (size_t i1 = 0; i1 < n1 - padX1; ++i1) {
    for (size_t i3 = 0; i3 < n3 - padY1; ++i3) {
      for (size_t i4 = 0; i4 < n4 - padY2; ++i4) {
        mat[i1 + 3][0][i3 + 3][i4 + 3] = mat[i1 + 3][1][i3 + 3][i4 + 3] =
            mat[i1 + 3][2][i3 + 3][i4 + 3] = mat[i1 + 3][3][i3 + 3][i4 + 3];
        mat[i1 + 3][n2 - 1][i3 + 3][i4 + 3] =
            mat[i1 + 3][n2 - 2][i3 + 3][i4 + 3] =
                mat[i1 + 3][n2 - 3][i3 + 3][i4 + 3] =
                    mat[i1 + 3][n2 - 4][i3 + 3][i4 + 3];
      }
    }
  }

  for (size_t i2 = 0; i2 < n2 - padX1; ++i2) {
    for (size_t i3 = 0; i3 < n3 - padY1; ++i3) {
      for (size_t i4 = 0; i4 < n4 - padY2; ++i4) {
        mat[0][i2 + 3][i3 + 3][i4 + 3] = mat[1][i2 + 3][i3 + 3][i4 + 3] =
            mat[2][i2 + 3][i3 + 3][i4 + 3] = mat[3][i2 + 3][i3 + 3][i4 + 3];
        mat[n1 - 1][i2 + 3][i3 + 3][i4 + 3] =
            mat[n1 - 2][i2 + 3][i3 + 3][i4 + 3] =
                mat[n1 - 3][i2 + 3][i3 + 3][i4 + 3] =
                    mat[n1 - 4][i2 + 3][i3 + 3][i4 + 3];
      }
    }
  }

  for (size_t i1 = 0; i1 < 3; ++i1) {
    for (size_t i2 = 0; i2 < 3; ++i2) {
      for (size_t i3 = 0; i3 < 3; ++i3) {
        for (size_t i4 = 0; i4 < 3; ++i4) {
          mat[i1][i2][i3][i4] = mat[3][3][3][3];

          mat[n1 - 3 + i1][i2][i3][i4] = mat[n1 - 4][3][3][3];

          mat[i1][n2 - 3 + i2][i3][i4] = mat[3][n2 - 4][3][3];
          mat[n1 - 3 + i1][n2 - 3 + i2][i3][i4] = mat[n1 - 4][n2 - 4][3][3];

          mat[i1][i2][n3 - 3 + i3][i4] = mat[3][3][n3 - 4][3];
          mat[n1 - 3 + i1][i2][n3 - 3 + i3][i4] = mat[n1 - 4][3][n3 - 4][3];
          mat[i1][n2 - 3 + i2][n3 - 3 + i3][i4] = mat[3][n2 - 4][n3 - 4][3];
          mat[n1 - 3 + i1][n2 - 3 + i2][n3 - 3 + i3][i4] =
              mat[n1 - 4][n2 - 4][n3 - 4][3];

          mat[i1][i2][i3][n4 - 3 + i4] = mat[3][3][3][n4 - 4];
          mat[n1 - 3 + i1][i2][i3][n4 - 3 + i4] = mat[n1 - 4][3][3][n4 - 4];
          mat[i1][n2 - 3 + i2][i3][n4 - 3 + i4] = mat[3][n2 - 4][3][n4 - 4];
          mat[n1 - 3 + i1][n2 - 3 + i2][i3][n4 - 3 + i4] =
              mat[n1 - 4][n2 - 4][3][n4 - 4];
          mat[i1][i2][n3 - 3 + i3][n4 - 3 + i4] = mat[3][3][n3 - 4][n4 - 4];
          mat[n1 - 3 + i1][i2][n3 - 3 + i3][n4 - 3 + i4] =
              mat[n1 - 4][3][n3 - 4][n4 - 4];
          mat[i1][n2 - 3 + i2][n3 - 3 + i3][n4 - 3 + i4] =
              mat[3][n2 - 4][n3 - 4][n4 - 4];
          mat[n1 - 3 + i1][n2 - 3 + i2][n3 - 3 + i3][n4 - 3 + i4] =
              mat[n1 - 4][n2 - 4][n3 - 4][n4 - 4];
        }
      }
    }
  }
}

tuple<double, double, double> stencils(const double &v1, const double &v2,
                                       const double &v3, const double &v4,
                                       const double &v5) noexcept {
  // smoothness of stencils
  auto S1 = 13. / 12 * (v1 - 2 * v2 + v3) * (v1 - 2 * v2 + v3) +
            1. / 4 * (v1 - 4 * v2 + 3 * v3) * (v1 - 4 * v2 + 3 * v3);
  auto S2 = 13. / 12 * (v2 - 2 * v3 + v4) * (v2 - 2 * v3 + v4) +
            1. / 4 * (v2 - v4) * (v2 - v4);
  auto S3 = 13. / 12 * (v3 - 2 * v4 + v5) * (v3 - 2 * v4 + v5) +
            1. / 4 * (3 * v3 - 4 * v4 + v5) * (3 * v3 - 4 * v4 + v5);
  // return make_tuple(S1, S2, S3);
  return {S1, S2, S3};
}

double derf(const double &v1, const double &v2, const double &v3,
            const double &v4, const double &v5) noexcept {
  // smoothness of stencils
  auto S1 = 13. / 12. * (v1 - 2. * v2 + v3) * (v1 - 2. * v2 + v3) +
            1. / 4. * (v1 - 4. * v2 + 3. * v3) * (v1 - 4. * v2 + 3. * v3);
  auto S2 = 13. / 12. * (v2 - 2. * v3 + v4) * (v2 - 2. * v3 + v4) +
            1. / 4. * (v2 - v4) * (v2 - v4);
  auto S3 = 13. / 12. * (v3 - 2. * v4 + v5) * (v3 - 2. * v4 + v5) +
            1. / 4. * (3. * v3 - 4. * v4 + v5) * (3. * v3 - 4. * v4 + v5);
  // find the maximum of vi's
  auto maxv2 = max({v1 * v1, v2 * v2, v3 * v3, v4 * v4, v5 * v5});
  auto eps1 = 1e-6 * maxv2 + 1e-99;
  // weights
  auto a1 = .1 / ((S1 + eps1) * (S1 + eps1));
  auto a2 = .6 / ((S2 + eps1) * (S2 + eps1));
  auto a3 = .3 / ((S3 + eps1) * (S3 + eps1));
  auto w1 = a1 / (a1 + a2 + a3);
  auto w2 = a2 / (a1 + a2 + a3);
  auto w3 = a3 / (a1 + a2 + a3);
  auto f1x = v1 / 3. - 7. * v2 / 6 + 11. * v3 / 6.;
  auto f2x = -v2 / 6. + 5. * v3 / 6. + v4 / 3.;
  auto f3x = v3 / 3. + 5. * v4 / 6. - v5 / 6.;
  return w1 * f1x + w2 * f2x + w3 * f3x;
}

template <size_t n1, size_t n2, size_t n3, size_t n4>
derivativeApproximations<n1 - 6, n2 - 6, n3 - 6, n4 - 6>
derf(const Array4D<n1, n2, n3, n4> &mat) noexcept {
  Array4D<n1 - 1, n2 - 6, n3 - 6, n4 - 6> wx1;
  Array4D<n1 - 6, n2 - 1, n3 - 6, n4 - 6> wx2;
  Array4D<n1 - 6, n2 - 6, n3 - 1, n4 - 6> wy1;
  Array4D<n1 - 6, n2 - 6, n3 - 6, n4 - 1> wy2;

  for (size_t i1 = 0; i1 < n1 - 1; ++i1) {
    for (size_t i2 = 0; i2 < n2 - 6; ++i2) {
      for (size_t i3 = 0; i3 < n3 - 6; ++i3) {
        for (size_t i4 = 0; i4 < n4 - 6; ++i4) {
          wx1[i1][i2][i3][i4] = (mat[i1 + 1][i2 + 3][i3 + 3][i4 + 3] -
                                 mat[i1][i2 + 3][i3 + 3][i4 + 3]) /
                                dx1;
        }
      }
    }
  }
  for (size_t i1 = 0; i1 < n1 - 6; ++i1) {
    for (size_t i2 = 0; i2 < n2 - 1; ++i2) {
      for (size_t i3 = 0; i3 < n3 - 6; ++i3) {
        for (size_t i4 = 0; i4 < n4 - 6; ++i4) {
          wx2[i1][i2][i3][i4] = (mat[i1 + 3][i2 + 1][i3 + 3][i4 + 3] -
                                 mat[i1 + 3][i2][i3 + 3][i4 + 3]) /
                                dx2;
        }
      }
    }
  }
  for (size_t i1 = 0; i1 < n1 - 6; ++i1) {
    for (size_t i2 = 0; i2 < n2 - 6; ++i2) {
      for (size_t i3 = 0; i3 < n3 - 1; ++i3) {
        for (size_t i4 = 0; i4 < n4 - 6; ++i4) {
          wy1[i1][i2][i3][i4] = (mat[i1 + 3][i2 + 3][i3 + 1][i4 + 3] -
                                 mat[i1 + 3][i2 + 3][i3][i4 + 3]) /
                                dy1;
        }
      }
    }
  }
  for (size_t i1 = 0; i1 < n1 - 6; ++i1) {
    for (size_t i2 = 0; i2 < n2 - 6; ++i2) {
      for (size_t i3 = 0; i3 < n3 - 6; ++i3) {
        for (size_t i4 = 0; i4 < n4 - 1; ++i4) {
          wy2[i1][i2][i3][i4] = (mat[i1 + 3][i2 + 3][i3 + 3][i4 + 1] -
                                 mat[i1 + 3][i2 + 3][i3 + 3][i4]) /
                                dy2;
        }
      }
    }
  }

  derivativeApproximations<n1 - 6, n2 - 6, n3 - 6, n4 - 6> result;

  for (size_t i1 = 0; i1 < n1 - 6; ++i1) {
    for (size_t i2 = 0; i2 < n2 - 6; ++i2) {
      for (size_t i3 = 0; i3 < n3 - 6; ++i3) {
        for (size_t i4 = 0; i4 < n4 - 6; ++i4) {
          result.fx1m[i1][i2][i3][i4] =
              derf(wx1[i1][i2][i3][i4], wx1[i1 + 1][i2][i3][i4],
                   wx1[i1 + 2][i2][i3][i4], wx1[i1 + 3][i2][i3][i4],
                   wx1[i1 + 4][i2][i3][i4]);
          result.fx1p[i1][i2][i3][i4] =
              derf(wx1[i1 + 5][i2][i3][i4], wx1[i1 + 4][i2][i3][i4],
                   wx1[i1 + 3][i2][i3][i4], wx1[i1 + 2][i2][i3][i4],
                   wx1[i1 + 1][i2][i3][i4]);

          result.fx2m[i1][i2][i3][i4] =
              derf(wx2[i1][i2][i3][i4], wx2[i1][i2 + 1][i3][i4],
                   wx2[i1][i2 + 2][i3][i4], wx2[i1][i2 + 3][i3][i4],
                   wx2[i1][i2 + 4][i3][i4]);
          result.fx2p[i1][i2][i3][i4] =
              derf(wx2[i1][i2 + 5][i3][i4], wx2[i1][i2 + 4][i3][i4],
                   wx2[i1][i2 + 3][i3][i4], wx2[i1][i2 + 2][i3][i4],
                   wx2[i1][i2 + 1][i3][i4]);

          result.fy1m[i1][i2][i3][i4] =
              derf(wy1[i1][i2][i3][i4], wy1[i1][i2][i3 + 1][i4],
                   wy1[i1][i2][i3 + 2][i4], wy1[i1][i2][i3 + 3][i4],
                   wy1[i1][i2][i3 + 4][i4]);
          result.fy1p[i1][i2][i3][i4] =
              derf(wy1[i1][i2][i3 + 5][i4], wy1[i1][i2][i3 + 4][i4],
                   wy1[i1][i2][i3 + 3][i4], wy1[i1][i2][i3 + 2][i4],
                   wy1[i1][i2][i3 + 1][i4]);

          result.fy2m[i1][i2][i3][i4] =
              derf(wy2[i1][i2][i3][i4], wy2[i1][i2][i3][i4 + 1],
                   wy2[i1][i2][i3][i4 + 2], wy2[i1][i2][i3][i4 + 3],
                   wy2[i1][i2][i3][i4 + 4]);
          result.fy2p[i1][i2][i3][i4] =
              derf(wy2[i1][i2][i3][i4 + 5], wy2[i1][i2][i3][i4 + 4],
                   wy2[i1][i2][i3][i4 + 3], wy2[i1][i2][i3][i4 + 2],
                   wy2[i1][i2][i3][i4 + 1]);
        }
      }
    }
  }

  return result;
}

double hamiltonianInDirection(const double &fp, const double &fm,
                              const double &changeInSign, const double &interm,
                              const double &a) noexcept {
  double res = 0;
  res += (1 - changeInSign) *
         ((interm < 1e-99) * interm * fp + (interm > 1e-99) * interm * fm);
  res += changeInSign * interm * (fp + fm) / 2.;
  res -= changeInSign * a * (fp - fm) / 2.;
  return res;
}

template <size_t n1, size_t n2, size_t n3, size_t n4>
Array4D<n1, n2, n3, n4>
hamiltonian(const derivativeApproximations<n1, n2, n3, n4> &approx) noexcept {
  Array4D<n1, n2, n3, n4> hamiltonian;

  for (size_t i1 = 0; i1 < n1; ++i1) {
    for (size_t i2 = 0; i2 < n2; ++i2) {
      for (size_t i3 = 0; i3 < n3; ++i3) {
        for (size_t i4 = 0; i4 < n4; ++i4) {
          const auto fx1p = approx.fx1p[i1][i2][i3][i4];
          const auto fx1m = approx.fx1m[i1][i2][i3][i4];
          const auto chooseAx1 = wwx1[i1][i2][i3][i4] * (fx1p * fx1m < 1e-99);
          const auto intermX1 =
              hx1[i1][i2][i3][i4] +
              (mx1[i1][i2][i3][i4] * (fx1p + fx1m) / 2. > 1e-99) *
                  mx1[i1][i2][i3][i4];
          hamiltonian[i1][i2][i3][i4] =
              (1 - chooseAx1) * ((intermX1 < 1e-99) * intermX1 * fx1p +
                                 (intermX1 > 1e-99) * intermX1 * fx1m);
          hamiltonian[i1][i2][i3][i4] +=
              chooseAx1 * intermX1 * (fx1p + fx1m) / 2.;
          hamiltonian[i1][i2][i3][i4] -=
              chooseAx1 * ax1[i1][i2][i3][i4] * (fx1p - fx1m) / 2.;

          const auto fx2p = approx.fx2p[i1][i2][i3][i4];
          const auto fx2m = approx.fx2m[i1][i2][i3][i4];
          const auto chooseAx2 = wwx2[i1][i2][i3][i4] * (fx2p * fx2m < 1e-99);
          const auto intermX2 =
              hx2[i1][i2][i3][i4] +
              (mx2[i1][i2][i3][i4] * (fx2p + fx2m) / 2. > 1e-99) *
                  mx2[i1][i2][i3][i4];
          hamiltonian[i1][i2][i3][i4] +=
              (1 - chooseAx2) * ((intermX2 < 1e-99) * intermX2 * fx2p +
                                 (intermX2 > 1e-99) * intermX2 * fx2m);
          hamiltonian[i1][i2][i3][i4] +=
              chooseAx2 * intermX2 * (fx2p + fx2m) / 2.;
          hamiltonian[i1][i2][i3][i4] -=
              chooseAx2 * ax2[i1][i2][i3][i4] * (fx2p - fx2m) / 2.;

          const auto fy1p = approx.fy1p[i1][i2][i3][i4];
          const auto fy1m = approx.fy1m[i1][i2][i3][i4];
          const auto chooseAy1 = wwy1[i1][i2][i3][i4] * (fy1p * fy1m < 1e-99);
          const auto intermY1 =
              hy1[i1][i2][i3][i4] +
              (my1[i1][i2][i3][i4] * (fy1p + fy1m) / 2. > 1e-99) *
                  my1[i1][i2][i3][i4];
          hamiltonian[i1][i2][i3][i4] +=
              (1 - chooseAy1) * ((intermY1 < 1e-99) * intermY1 * fy1p +
                                 (intermY1 > 1e-99) * intermY1 * fy1m);
          hamiltonian[i1][i2][i3][i4] +=
              chooseAy1 * intermY1 * (fy1p + fy1m) / 2.;
          hamiltonian[i1][i2][i3][i4] -=
              chooseAy1 * ay1[i1][i2][i3][i4] * (fy1p - fy1m) / 2.;

          const auto fy2p = approx.fy2p[i1][i2][i3][i4];
          const auto fy2m = approx.fy2m[i1][i2][i3][i4];
          const auto chooseAy2 = wwy2[i1][i2][i3][i4] * (fy2p * fy2m < 1e-99);
          const auto intermY2 =
              hy2[i1][i2][i3][i4] +
              (my2[i1][i2][i3][i4] * (fy2p + fy2m) / 2. > 1e-99) *
                  my2[i1][i2][i3][i4];
          hamiltonian[i1][i2][i3][i4] +=
              (1 - chooseAy2) * ((intermY2 < 1e-99) * intermY2 * fy2p +
                                 (intermY2 > 1e-99) * intermY2 * fy2m);
          hamiltonian[i1][i2][i3][i4] +=
              chooseAy2 * intermY2 * (fy2p + fy2m) / 2.;
          hamiltonian[i1][i2][i3][i4] -=
              chooseAy2 * ay2[i1][i2][i3][i4] * (fy2p - fy2m) / 2.;
        }
      }
    }
  }

  return hamiltonian;
}

template <size_t n1, size_t n2, size_t n3, size_t n4>
void performTimeStep(Array4D<n1, n2, n3, n4> &tp1Matrix,
                     const Array4D<n1, n2, n3, n4> &tMatrix) {
  auto approx = derf(tMatrix);
  auto hamilton = hamiltonian(approx);
  for (size_t i1 = 0; i1 < n1 - padX1; ++i1) {
    for (size_t i2 = 0; i2 < n2 - padX2; ++i2) {
      for (size_t j1 = 0; j1 < n3 - padY1; ++j1) {
        for (size_t j2 = 0; j2 < n4 - padY2; ++j2) {
          tp1Matrix[i1 + 3][i2 + 3][j1 + 3][j2 + 3] =
              max((1 - L * dt) * tMatrix[i1 + 3][i2 + 3][j1 + 3][j2 + 3] -
                      hamilton[i1][i2][j1][j2] * dt,
                  distSgnd[i1][i2][j1][j2]);
        }
      }
    }
  }
  matrixContinuation(tp1Matrix);
}

int main(int argc, char *argv[]) {
  // generate(x1.begin(), x1.end(),
  //          [n = 0]() mutable { return (-padX1 + n++) * dx1; });
  // x1[padX1] = 0;
  // x1[padX1 + Nx1 - 1] = Imax1;
  // generate(y1.begin(), y1.end(),
  //          [n = 0]() mutable { return (-padY1 + n++) * dy1; });
  // y1[padY1] = 0;
  // y1[padY1 + Ny1 - 1] = ymax1;
  // generate(x2.begin(), x2.end(),
  //          [n = 0]() mutable { return (-padX2 + n++) * dx2; });
  // x2[padX2] = 0;
  // x2[padX2 + Nx2 - 1] = Imax2;
  // generate(y2.begin(), y2.end(),
  //          [n = 0]() mutable { return (-padY2 + n++) * dy2; });
  // y2[padY2] = 0;
  // y2[padY2 + Ny2 - 1] = ymax2;

  generate(x1.begin(), x1.end(), [n = 0]() mutable { return (n++) * dx1; });
  x1[Nx1 - 1] = Imax1;
  generate(x2.begin(), x2.end(), [n = 0]() mutable { return (n++) * dx2; });
  x2[Nx2 - 1] = Imax2;
  generate(y1.begin(), y1.end(), [n = 0]() mutable { return (n++) * dy1; });
  y1[Ny1 - 1] = ymax1;
  generate(y2.begin(), y2.end(), [n = 0]() mutable { return (n++) * dy2; });
  y2[Ny2 - 1] = ymax2;

  // compute the partial derivatives
  for (size_t i1 = 0; i1 < x1.size(); ++i1) {
    for (size_t i2 = 0; i2 < x2.size(); ++i2) {
      for (size_t j1 = 0; j1 < y1.size(); ++j1) {
        for (size_t j2 = 0; j2 < y2.size(); ++j2) {
          hx1[i1][i2][j1][j2] =
              Gamma1 * x1[i1] - (1 - x1[i1]) * (b11 * y1[j1] + b12 * y2[j2]);
          hx2[i1][i2][j1][j2] =
              Gamma2 * x2[i2] - (1 - x2[i2]) * (b21 * y1[j1] + b22 * y2[j2]);
          hy1[i1][i2][j1][j2] =
              mu1 * y1[j1] - (1 - y1[j1]) * (c11 * x1[i1] + c12 * x2[i2]);
          hy2[i1][i2][j1][j2] =
              mu2 * y2[j2] - (1 - y2[j2]) * (c21 * x1[i1] + c22 * x2[i2]);

          mx1[i1][i2][j1][j2] =
              k * umax1 * (1 - x1[i1]) * (b11 * y1[j1] + b12 * y2[j2]);
          mx2[i1][i2][j1][j2] =
              k * umax2 * (1 - x2[i2]) * (b21 * y1[j1] + b22 * y2[j2]);
          my1[i1][i2][j1][j2] =
              k * (1 - y1[j1]) * (c11 * umax1 * x1[i1] + c12 * umax2 * x2[i2]);
          my2[i1][i2][j1][j2] =
              k * (1 - y2[j2]) * (c21 * umax1 * x1[i1] + c22 * umax2 * x2[i2]);

          // find maximum range of gridpoints
          wwx1[i1][i2][j1][j2] =
              (hx1[i1][i2][j1][j2] *
               (hx1[i1][i2][j1][j2] + mx1[i1][i2][j1][j2])) < 1e-99;
          ax1[i1][i2][j1][j2] =
              max(abs(hx1[i1][i2][j1][j2]),
                  abs(hx1[i1][i2][j1][j2] + mx1[i1][i2][j1][j2]));
          wwx2[i1][i2][j1][j2] =
              (hx1[i1][i2][j1][j2] *
               (hx2[i1][i2][j1][j2] + mx2[i1][i2][j1][j2])) < 1e-99;
          ax2[i1][i2][j1][j2] =
              max(abs(hx2[i1][i2][j1][j2]),
                  abs(hx2[i1][i2][j1][j2] + mx2[i1][i2][j1][j2]));
          wwy1[i1][i2][j1][j2] =
              (hy1[i1][i2][j1][j2] *
               (hy1[i1][i2][j1][j2] + my1[i1][i2][j1][j2])) < 1e-99;
          ay1[i1][i2][j1][j2] =
              max(abs(hy1[i1][i2][j1][j2]),
                  abs(hy1[i1][i2][j1][j2] + my1[i1][i2][j1][j2]));
          wwy2[i1][i2][j1][j2] =
              (hy1[i1][i2][j1][j2] *
               (hy1[i1][i2][j1][j2] + my1[i1][i2][j1][j2])) < 1e-99;
          ay2[i1][i2][j1][j2] =
              max(abs(hy2[i1][i2][j1][j2]),
                  abs(hy2[i1][i2][j1][j2] + my2[i1][i2][j1][j2]));

          // set signed distance
          distSgnd[i1][i2][j1][j2] =
              -min(abs(x1[i1] - Imax1), abs(x2[i2] - Imax2));
          // -__builtin_sqrt((x1[i1] - Imax1) * (x1[i1] - Imax1) +
          //                 (x2[i2] - Imax2) * (x2[i2] - Imax2));

          // set initial approximation
          w[i1 + 3][i2 + 3][j1 + 3][j2 + 3] =
              (__builtin_sin((x1[i1] + x2[i2] - .5) * M_PI * 2.)) * 1e-6 +
              distSgnd[i1][i2][j1][j2];
        }
      }
    }
  }

  matrixContinuation(w);
  auto wnew0 = w;
  auto wnew1 = w;

  auto currentError = 1.;
  vector<double> err;
  auto ncounter = 0;
  auto start = chrono::high_resolution_clock::now();

  while (currentError > errtol) {

    // Heun's predictor-corrector method
    performTimeStep(wnew0, w);
    performTimeStep(wnew1, wnew0);
    const auto wnew = halfSum(w, wnew1);

    currentError = 0;
    for (size_t i1 = 0; i1 < x1.size(); ++i1) {
      for (size_t i2 = 0; i2 < x2.size(); ++i2) {
        for (size_t j1 = 0; j1 < y1.size(); ++j1) {
          for (size_t j2 = 0; j2 < y2.size(); ++j2) {
            currentError =
                max(currentError, abs(w[i1 + 3][i2 + 3][j1 + 3][j2 + 3] -
                                      wnew[i1 + 3][i2 + 3][j1 + 3][j2 + 3]));
          }
        }
      }
    }
    err.push_back(currentError);
    if (!err.empty() && currentError > err.back()) {
      cout << "AAAA" << "\n";
      cout << "Steps: " << ncounter << "\n";
      cout << "Error last: " << err.back() << "Current error: " << currentError
           << "\n";
    }
    w = wnew;
    ncounter++;
    if (ncounter % 100 == 0) {
      cout << "Steps: " << ncounter << "\n";
      cout << "Error: " << err.back() << "\n";
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
  printSolution(w, finalSolution);

  // system(format("cat {}", outputFile.c_str()).c_str());
  // ifstream finalSolutionReRead(outputFile);
  // cout << string((istreambuf_iterator<char>(finalSolutionReRead)),
  //                istreambuf_iterator<char>());
}
