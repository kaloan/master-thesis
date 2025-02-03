#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <iterator>
#include <tuple>
// #include <numeric>
// #include <functional>
using namespace std;

template <size_t rows, size_t cols> struct derivativeApproximations {
  array<array<double, cols - 6>, rows - 6> fxp;
  array<array<double, cols - 6>, rows - 6> fxm;
  array<array<double, cols - 6>, rows - 6> fyp;
  array<array<double, cols - 6>, rows - 6> fym;
};

// parameters from Goldberg, Agusto 2021
// human
constexpr auto Gamma = 1. / 14;
// constexpr auto alpha = 0.082*0.5*exp(-1/3); // Botsuana
const auto alpha = 0.241 * 0.5 * exp(-1 / 3); // Zimbabwe
// mosquito
// beta = 0.082*0.1; mu=1/30; // Botsuana
constexpr auto Beta = 0.241 * 0.1; // Zimbabwe
constexpr auto mu = 1. / 10;       // Zimbabwe

// ratio mosquito : human
constexpr auto nu = 10.;
// maximum infected
constexpr auto Imax = .1;
// set to reduce the size of the domain where the viability kernel lies
constexpr auto ymax = .2;
// control
constexpr auto k = .6;
constexpr auto umax = .6;

// threshold for I at maximum control
const auto Is =
    (alpha * nu * (1 - k * umax) * (1 - k * umax) * Beta - Gamma * mu) /
    ((Gamma + alpha * nu * (1 - k * umax)) * (1 - k * umax) * Beta);

// predefine
constexpr auto ntim = 650.;
constexpr auto errtol = 1e-10;

// lambda - lipschitz constant for the system
const auto ll1 = 2 * (alpha * alpha * nu * nu + Beta * Beta + mu * mu);
const auto ll2 = sqrt(ll1);
const auto ll = 1.05 * ll2;

// maximum of partial derivatives, Hamiltonian
const array<double, 2> mpd = {1 + Gamma + alpha * nu * (1 + k * umax),
                              1 + mu + Beta * (1 + k * umax)};
// mpd = mpd + 1; ????

// numerical space step
constexpr auto Nx = 41;
constexpr auto dx = Imax / Nx;
constexpr auto Ny = 101;
constexpr auto dy = ymax / Ny;

// numerical time step dt < 1/ll (!)
// numerical time step
// const auto dt = 0.9 * min({1. / (mpd[0] / dx + mpd[1] / dy), 1. / ll});
const auto dt = 0.9 / max(mpd[0] / dx + mpd[1] / dy, ll);

template <typename Type, std::size_t... sizes>
auto concatenate(const std::array<Type, sizes> &...arrays) {
  std::array<Type, (sizes + ...)> result;
  std::size_t index{};

  ((std::copy_n(arrays.begin(), sizes, result.begin() + index), index += sizes),
   ...);

  return result;
}

template <typename T, size_t rows, size_t cols>
array<array<T, rows>, cols> transpose(const array<array<T, cols>, rows> &mat) {
  array<array<T, rows>, cols> transposed;

  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) {
      transposed[j][i] = mat[i][j];
    }
  }

  return transposed;
}

// x and y
constexpr auto padX = 6;
array<double, Nx + 2 * padX> x;
constexpr auto padY = 6;
array<double, Ny + 2 * padY> y;
// mdspan<double, extents<std::size_t, 1>> MDSPAN(x.data(), Nx + 2 * padX);
array<array<double, x.size()>, y.size()> h1;
array<array<double, x.size()>, y.size()> h2;
array<array<double, x.size()>, y.size()> m1;
array<array<double, x.size()>, y.size()> m2;
array<array<double, x.size()>, y.size()> wwx;
array<array<double, x.size()>, y.size()> wwy;
array<array<double, x.size()>, y.size()> ax;
array<array<double, x.size()>, y.size()> ay;
array<array<double, x.size()>, y.size()> gs;
// array<array<double, x.size()>, y.size()> wm;
array<array<double, x.size() + 6>, y.size() + 6> w;

template <typename T, size_t size>
array<T, size> sum(const array<T, size> &arr1, const array<T, size> &arr2) {
  array<T, size> res;
  for (size_t i = 0; i < size; i++) {
    res[i] = arr1[i] + arr2[i];
  }
  return res;
}

template <typename T, size_t size>
array<T, size> halfSum(const array<T, size> &arr1, const array<T, size> &arr2) {
  array<T, size> res;
  for (size_t i = 0; i < size; i++) {
    res[i] = (arr1[i] + arr2[i]) / 2.;
  }
  return res;
}

template <typename T, size_t rows, size_t cols>
array<array<T, cols>, rows> halfSum(const array<array<T, cols>, rows> &mat1,
                                    const array<array<T, cols>, rows> &mat2) {
  array<array<T, cols>, rows> res;
  for (size_t i = 0; i < rows; i++) {
    res[i] = halfSum(mat1[i], mat2[i]);
  }
  return res;
}

template <size_t cols, size_t rows>
double find_maxv2(const array<array<double, cols>, rows> &matrix) {
  auto maxv2 = 0.;
  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) {
      maxv2 = max(maxv2, matrix[i][j] * matrix[i][j]);
    }
  }
  return maxv2;
}

template <size_t cols, size_t rows>
array<array<double, cols>, rows>
max(const array<array<double, cols>, rows> &mat1,
    const array<array<double, cols>, rows> &mat2) {
  array<array<double, cols>, rows> result;
  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < cols; j++) {
      result[i][j] = max(mat1[i][j], mat2[i][j]);
    }
  }
  return result;
}

template <typename T, size_t size> void printArray(const array<T, size> &arr) {
  std::copy(std::begin(arr), std::end(arr),
            std::ostream_iterator<T>(std::cout, ", "));
  std::cout << std::endl;
}

template <typename T, size_t rows, size_t cols>
void printMatrix(const array<array<T, cols>, rows> &arr) {
  for (auto row : arr) {
    printArray(row);
  }
}

template <typename T, size_t rows, size_t cols>
void matrixContinuation(array<array<T, cols>, rows> &mat) {
  for (size_t j = 0; j < x.size(); j++) {
    mat[0][j + 3] = mat[1][j + 3] = mat[2][j + 3] = mat[3][j + 3];
    mat[y.size() + 3][j + 3] = mat[y.size() + 4][j + 3] =
        mat[y.size() + 5][j + 3] = mat[y.size() + 2][j + 3];
  }
  for (size_t i = 0; i < y.size(); i++) {
    mat[i + 3][0] = mat[i + 3][1] = mat[i + 3][2] = mat[i + 3][3];
    mat[i + 3][x.size() + 3] = mat[i + 3][x.size() + 4] =
        mat[i + 3][x.size() + 5] = mat[i + 3][x.size() + 2];
  }
}

tuple<double, double, double> stencils(const double &v1, const double &v2,
                                       const double &v3, const double &v4,
                                       const double &v5) {
  auto f1x = v1 / 3 - 7 * v2 / 6 + 11 * v3 / 6;
  auto f2x = -v2 / 6 + 5 * v3 / 6 + v4 / 3;
  auto f3x = v3 / 3 + 5 * v4 / 6 - v5 / 6;
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
            const double &v4, const double &v5, const double &eps1) {
  auto f1x = v1 / 3 - 7 * v2 / 6 + 11 * v3 / 6;
  auto f2x = -v2 / 6 + 5 * v3 / 6 + v4 / 3;
  auto f3x = v3 / 3 + 5 * v4 / 6 - v5 / 6;
  // smoothness of stencils
  auto S1 = 13. / 12 * (v1 - 2 * v2 + v3) * (v1 - 2 * v2 + v3) +
            1. / 4 * (v1 - 4 * v2 + 3 * v3) * (v1 - 4 * v2 + 3 * v3);
  auto S2 = 13. / 12 * (v2 - 2 * v3 + v4) * (v2 - 2 * v3 + v4) +
            1. / 4 * (v2 - v4) * (v2 - v4);
  auto S3 = 13. / 12 * (v3 - 2 * v4 + v5) * (v3 - 2 * v4 + v5) +
            1. / 4 * (3 * v3 - 4 * v4 + v5) * (3 * v3 - 4 * v4 + v5);
  // find the maximum of vi's element-wise
  // auto findMaxList = {v1 * v1, v2 * v2, v3 * v3, v4 * v4, v5 * v5};
  // auto maxv2 = *max_element(findMaxList.begin(), findMaxList.end());
  // auto eps1 = 1e-6 * maxv2 + 1e-99;
  // weights
  auto a1 = .1 / ((S1 + eps1) * (S1 + eps1));
  auto a2 = .6 / ((S2 + eps1) * (S2 + eps1));
  auto a3 = .3 / ((S3 + eps1) * (S3 + eps1));
  auto w1 = a1 / (a1 + a2 + a3);
  auto w2 = a2 / (a1 + a2 + a3);
  auto w3 = a3 / (a1 + a2 + a3);
  return w1 * f1x + w2 * f2x + w3 * f3x;
}

template <size_t n>
array<double, n> derf(const array<double, n> &v1, const array<double, n> &v2,
                      const array<double, n> &v3, const array<double, n> &v4,
                      const array<double, n> &v5, const double &eps1) {
  array<double, n> res;
  for (size_t i = 0; i < n; i++) {
    res[i] = derf(v1[i], v2[i], v3[i], v4[i], v5[i], eps1);
  }
  return res;
}

template <typename T, size_t rows, size_t cols>
derivativeApproximations<rows, cols>
derf(const array<array<T, cols>, rows> &mat) {
  array<array<double, cols - 6>, rows - 6> wx;
  array<array<double, cols - 6>, rows - 6> wy;
  for (int i = 0; i < rows - 6; i++) {
    for (int j = 0; j < cols - 6; j++) {
      wx[i][j] = (mat[i + 3][j + 4] - mat[i + 3][j + 3]) / dx;
      wy[i][j] = (mat[i + 4][j + 3] - mat[i + 3][j + 3]) / dy;
    }
  }

  derivativeApproximations<rows, cols> result;

  auto epsY = 1e-6 * find_maxv2(wy) + 1e-99;
  for (size_t i = 0; i < wy.size() - 5; i++) {
    result.fym[i] =
        derf(wy[i], wy[i + 1], wy[i + 2], wy[i + 3], wy[i + 4], epsY);
  }
  for (size_t i = 5; i < wy.size(); i++) {
    result.fyp[i - 4] =
        derf(wy[i], wy[i - 1], wy[i - 2], wy[i - 3], wy[i - 4], epsY);
  }

  auto epsX = 1e-6 * find_maxv2(wx) + 1e-99;
  auto transposedWx = transpose(wx);
  array<array<double, rows - 6>, cols - 6> transposedFx;
  for (int i = 0; i < cols - 6; i++) {
    for (int j = 0; j < rows - 6; j++) {
      transposedFx[i][j] = 0;
    }
  }
  for (size_t i = 0; i < transposedWx.size() - 5; i++) {
    transposedFx[i] =
        derf(transposedWx[i], transposedWx[i + 1], transposedWx[i + 2],
             transposedWx[i + 3], transposedWx[i + 4], epsX);
  }
  result.fxm = transpose(transposedFx);

  for (int i = 0; i < cols - 6; i++) {
    for (int j = 0; j < rows - 6; j++) {
      transposedFx[i][j] = 0;
    }
  }
  for (size_t i = 5; i < transposedWx.size(); i++) {
    transposedFx[i - 4] =
        derf(transposedWx[i], transposedWx[i - 1], transposedWx[i - 2],
             transposedWx[i - 3], transposedWx[i - 4], epsX);
  }
  result.fxp = transpose(transposedFx);

  return result;
}

template <size_t rows, size_t cols>
array<array<double, cols>, rows>
hamiltonian(const array<array<double, cols>, rows> &xp,
            const array<array<double, cols>, rows> &xm,
            const array<array<double, cols>, rows> &yp,
            const array<array<double, cols>, rows> &ym) {
  array<array<double, cols>, rows> hamiltonian;

  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      auto chooseAx = wwx[i][j] * (xp[i][j] * xm[i][j] < 1e-99);
      auto interm =
          h1[i][j] + (m1[i][j] * (xp[i][j] + xm[i][j]) / 2 > 1e-99) * m1[i][j];
      // upwind outside ww1 Osher & Shu (2.11) - (2.11) 1st case
      hamiltonian[i][j] +=
          (1 - chooseAx) * ((interm < 1e-99) * interm * xp[i][j] +
                            (interm > 1e-99) * interm * xm[i][j]);
      hamiltonian[i][j] += chooseAx * interm * (xp[i][j] + xm[i][j]) /
                           2; // Osher & Shu (2.11) -2nd case
      hamiltonian[i][j] -= chooseAx * ax[i][j] * (xp[i][j] - xm[i][j]) / 2;

      auto chooseAy = wwy[i][j] * (yp[i][j] * ym[i][j] < 1e-99);
      auto interm2 =
          h2[i][j] + (m2[i][j] * (yp[i][j] + ym[i][j]) / 2 > 1e-99) * m2[i][j];
      // upwind outside ww1 Osher & Shu (2.11) - (2.11) 1st case
      hamiltonian[i][j] +=
          (1 - chooseAy) * ((interm2 < 1e-99) * interm2 * yp[i][j] +
                            (interm2 > 1e-99) * interm2 * ym[i][j]);
      hamiltonian[i][j] += chooseAy * interm2 * (yp[i][j] + ym[i][j]) /
                           2; // Osher & Shu (2.11) -2nd case
      hamiltonian[i][j] -= chooseAy * ay[i][j] * (yp[i][j] - ym[i][j]) / 2;
    }
  }

  return hamiltonian;
}

int main(int argc, char *argv[]) {
  generate(x.begin(), x.end(),
           [n = 0]() mutable { return -(padX * dx) + (n++ * dx); });
  x[padX] = 0;
  x[Nx] = Imax;
  generate(y.begin(), y.end(),
           [n = 0]() mutable { return -(padY * dy) + (n++ * dy); });
  y[padY] = 0;
  y[Ny] = ymax;
  cout << "XXXXXX:";
  printArray(x);
  cout << "YYYYYY:";
  printArray(y);
  // compute the partial derivatives
  for (size_t i = 0; i < y.size(); i++) {
    for (size_t j = 0; j < x.size(); j++) {
      h1[i][j] = Gamma * x[j] - alpha * nu * (1 - x[j]) * y[i];
      h2[i][j] = mu * y[i] - Beta * x[j] * (1 - y[i]);
      m1[i][j] = alpha * nu * k * umax * (1 - x[j]) * y[i];
      m2[i][j] = Beta * k * umax * x[j] * (1 - y[i]);

      // find maximum range of gridpoints
      // where a_x must be chosen according to LLF
      wwx[i][j] = (h1[i][j] * (h1[i][j] + m1[i][j])) < 1e-99;
      ax[i][j] = max(abs(h1[i][j]), abs(h1[i][j] + m1[i][j]));

      // find maximum range of gridpoints
      // where a_y must be chosen according to LLF
      wwy[i][j] = (h2[i][j] * (h2[i][j] + m2[i][j])) < 1e-99;
      ay[i][j] = max(abs(h2[i][j]), abs(h2[i][j] + m2[i][j]));

      gs[i][j] = x[j] - Imax;
      // Initial guess
      w[i + 3][j + 3] = (sin((x[j] - .5) * M_PI * 2)) * 1e-6 + gs[i][j];
    }
  }
  matrixContinuation(w);
  auto wnew0 = w;
  auto wnew1 = w;
  // printMatrix(w);
  auto err0 = 1.;
  while (err0 > errtol) {
    // Heun's predictor-corrector method
    auto firstApprox = derf(w);
    auto hamilton = hamiltonian(firstApprox.fxm, firstApprox.fxp,
                                firstApprox.fym, firstApprox.fyp);
    for (size_t i = 0; i < y.size(); i++) {
      for (size_t j = 0; j < x.size(); j++) {
        wnew0[i + 3][j + 3] = max(
            (1 - ll * dt) * w[i + 3][j + 3] - hamilton[i][j] * dt, gs[i][j]);
      }
    }
    matrixContinuation(wnew0);

    auto secondApprox = derf(wnew0);
    hamilton = hamiltonian(secondApprox.fxm, secondApprox.fxp, secondApprox.fym,
                           secondApprox.fyp);
    for (size_t i = 0; i < y.size(); i++) {
      for (size_t j = 0; j < x.size(); j++) {
        wnew1[i + 3][j + 3] =
            max((1 - ll * dt) * wnew0[i + 3][j + 3] - hamilton[i][j] * dt,
                gs[i][j]);
      }
    }
    matrixContinuation(wnew1);

    auto wnew = halfSum(wnew0, wnew1);
    err0 = 0;
    for (size_t i = 0; i < y.size(); i++) {
      for (size_t j = 0; j < x.size(); j++) {
        err0 = max(err0, abs(w[i + 3][j + 3] - wnew[i + 3][j + 3]));
      }
    }
    w = wnew;
  }
}
