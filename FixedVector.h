#pragma once

#include <array>

template <size_t N>
requires (N > 0)
struct FixedVector {

public:
  std::array<double, N> coords;

  inline FixedVector() noexcept : coords{} {}

  inline FixedVector(const std::array<double, N> &coords) noexcept : coords(coords) {}

  inline FixedVector(std::array<double, N> &&coords) noexcept : coords(coords) {}

  inline constexpr double operator[](size_t index) noexcept {
    return coords[index];
  }

};
