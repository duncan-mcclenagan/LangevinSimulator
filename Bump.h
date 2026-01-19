#pragma once

#include <cmath>

struct Bump {
  const double center;
  const double width;
  const double height;

  constexpr Bump() : center(0), width(0), height(0) {}

  constexpr Bump(double center, double width, double height)
    : center(center), width(width > 0 ? width : -width), height(height > 0 ? height : -height) {
  }

  inline bool inBump(double position) const {
    return std::abs(position - center) < width;
  }
};