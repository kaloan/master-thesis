#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <iterator>
#include <tuple>
// #include <numeric>
// #include <functional>
using namespace std;

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
const array<double, 2> mpd = {Gamma + alpha * nu * (1 + k * umax),
                              mu + Beta * (1 + k * umax)};
// mpd = mpd + 1; ????

// numerical space step
constexpr auto Nx = 41;
constexpr auto dx = Imax / Nx;
constexpr auto Ny = 101;
constexpr auto dy = ymax / Ny;

// numerical time step dt < 1/ll (!)
// numerical time step
// const auto dt = 0.9 * min({1. / (mpd[0] / dx + mpd[1] / dy), 1. / ll});
const auto dt = 0.9 / max({mpd[0] / dx + mpd[1] / dy, ll});

template <typename Type, std::size_t... sizes>
auto concatenate(const std::array<Type, sizes> &...arrays) {
  std::array<Type, (sizes + ...)> result;
  std::size_t index{};

  ((std::copy_n(arrays.begin(), sizes, result.begin() + index), index += sizes),
   ...);

  return result;
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
void matrixContinuation(const array<array<T, cols>, rows> &mat) {
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
  auto S1 = 13. / 12 * pow((v1 - 2 * v2 + v3), 2) +
            1. / 4 * pow((v1 - 4 * v2 + 3 * v3), 2);
  auto S2 = 13. / 12 * pow((v2 - 2 * v3 + v4), 2) + 1. / 4 * pow((v2 - v4), 2);
  auto S3 = 13. / 12 * pow((v3 - 2 * v4 + v5), 2) +
            1. / 4 * pow((3 * v3 - 4 * v4 + v5), 2);
  // return make_tuple(S1, S2, S3);
  return {S1, S2, S3};
}

double derf(const double &v1, const double &v2, const double &v3,
            const double &v4, const double &v5, const double &eps1) {
  auto f1x = v1 / 3 - 7 * v2 / 6 + 11 * v3 / 6;
  auto f2x = -v2 / 6 + 5 * v3 / 6 + v4 / 3;
  auto f3x = v3 / 3 + 5 * v4 / 6 - v5 / 6;
  // smoothness of stencils
  auto S1 = 13. / 12 * pow((v1 - 2 * v2 + v3), 2) +
            1. / 4 * pow((v1 - 4 * v2 + 3 * v3), 2);
  auto S2 = 13. / 12 * pow((v2 - 2 * v3 + v4), 2) + 1. / 4 * pow((v2 - v4), 2);
  auto S3 = 13. / 12 * pow((v3 - 2 * v4 + v5), 2) +
            1. / 4 * pow((3 * v3 - 4 * v4 + v5), 2);
  // find the maximum of vi's element-wise
  // auto findMaxList = {v1 * v1, v2 * v2, v3 * v3, v4 * v4, v5 * v5};
  // auto maxv2 = *max_element(findMaxList.begin(), findMaxList.end());
  // auto eps1 = 1e-6 * maxv2 + 1e-99;
  // weights
  auto a1 = .1 / pow((S1 + eps1), 2);
  auto a2 = .6 / pow((S2 + eps1), 2);
  auto a3 = .3 / pow((S3 + eps1), 2);
  auto w1 = a1 / (a1 + a2 + a3);
  auto w2 = a2 / (a1 + a2 + a3);
  auto w3 = a3 / (a1 + a2 + a3);
  return w1 * f1x + w2 * f2x + w3 * f3x;
}

template <size_t n>
array<double, n> derf(array<double, n> v1, array<double, n> v2,
                      array<double, n> v3, array<double, n> v4,
                      array<double, n> v5) {
  auto maxv2 = 0.;
  for (size_t i = 0; i < n; i++) {
    auto findMaxList = {v1[i] * v1[i], v2[i] * v2[i], v3[i] * v3[i],
                        v4[i] * v4[i], v5[i] * v5[i]};
    auto maxv2Possible = *max_element(findMaxList.begin(), findMaxList.end());
    maxv2 = max(maxv2, maxv2Possible);
  }
  auto eps1 = 1e-6 * maxv2 + 1e-99;
  array<double, n> res;
  for (size_t i = 0; i < n; i++) {
    res[i] = derf(v1[i], v2[i], v3[i], v4[i], v5[i], eps1);
  }
  return res;
}

template <size_t rows, size_t cols> array<array<double, cols>, rows> a() {
  return {};
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
  // printMatrix(w);
  auto err0 = 1.;
  while (err0 > errtol) {
    // Heun's predictor-corrector method
    array<array<double, x.size() + 6>, y.size()> wy;
    array<array<double, x.size()>, y.size() + 6> wx;
    for (int i = 0; i < y.size(); i++) {
      for (int j = 0; j < x.size(); j++) {
        wy[i][j + 3] = w[i + 1][j + 3] - w[i][j + 3];
        wx[i + 3][j] = w[i + 3][j + 1] - w[i + 3][j];
      }
    }
  }
}
