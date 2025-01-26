#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <iterator>
#include <numeric>
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
// array<double, padX> xm;
// array<double, padX> xp;
// array<double, Nx> xb;
array<double, Nx + 2 * padX> x;
constexpr auto padY = 6;
array<double, Ny + 2 * padY> y;

template <typename T, size_t size>
array<T, size> sum(array<T, size> arr1, array<T, size> arr2) {
  array<T, size> res;
  for (size_t i = 0; i < size; i++) {
    res[i] = arr1[i] + arr2[i];
  }
  return res;
}

template <typename T, size_t size> void printArray(array<T, size> arr) {
  std::copy(std::begin(arr), std::end(arr),
            std::ostream_iterator<T>(std::cout, ", "));
  std::cout << std::endl;
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
}
