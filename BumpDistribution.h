#pragma once

#include <array>
#include <algorithm>
#include <cmath>
#include <optional>
#include <stdexcept>
#include "Bump.h"

template<int NBumps> requires (NBumps > 0)
class BumpDistribution {
private:
  const std::array<Bump, NBumps> bumps;
  const double inverseSegmentSpacing;
  static constexpr int centerSegmentIndex{ static_cast<int>(NBumps) / 2 };

public:
  constexpr BumpDistribution() : bumps(std::array<Bump, NBumps>()), inverseSegmentSpacing(0) {}

  constexpr BumpDistribution(std::array<Bump, NBumps> bumps, double segmentSpacing)
    : bumps(bumps), inverseSegmentSpacing(1 / (segmentSpacing > 0 ? segmentSpacing : -segmentSpacing)) {
  }

  Bump getAtSegmentIndex(int index) const {
    return bumps[index];
  }

  int getNearestSegmentIndex(double position) const {
    int index{ static_cast<int>(std::lround(position * inverseSegmentSpacing)) + centerSegmentIndex };
    return std::min(std::max(0, index), NBumps - 1);
  }

  Bump getNearestBump(double position) const {
    return bumps[getNearestSegmentIndex(position)];
  }

};
