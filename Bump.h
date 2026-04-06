#pragma once

#include <cmath>
#include "CoordinateSpan.h"

struct BumpData {
  const double* center;
  const double squareWidth;
  const double height;

  constexpr BumpData(const double* center, double width, double height)
    : center(center), squareWidth(squareWidth > 0 ? squareWidth : -squareWidth), height(height > 0 ? height : -height) {
  }

};

class Bump {
public:
  const CoordinateSpan<const double> center;
  const double squareWidth;
  const double height;

  constexpr Bump(CoordinateSpan<const double> center, double width, double height)
    : center(center), squareWidth(squareWidth > 0 ? squareWidth : -squareWidth), height(height > 0 ? height : -height) {
  }

  constexpr const void* uniqueIdentifier() const {
    return &center[0, 0];
  }

};
