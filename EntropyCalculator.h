#pragma once

#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include "PhaseCloud.h"
#include "BumpDistribution.h"
#include "Bump.h"

class EntropyCalculator {
private:
  std::unordered_map<int, int> bumpAccumulants{};

public:
  template<int NBumps>
  double calculateEntropy(const PhaseCloud& phaseCloud, const BumpDistribution<NBumps>& bumps) {
    if (phaseCloud.numberOfParticles != 1) {
      throw std::invalid_argument("Phase cloud number of particles must be one.");
    }
    if (phaseCloud.spatialDimension != 1) {
      throw std::invalid_argument("Phase cloud spatial dimension must be one.");
    }

    bumpAccumulants.clear();

    const double* position = phaseCloud.rawPositionConst();
    int totalInBumps{ 0 };
    for (int i{ 0 }; i < phaseCloud.numberOfPoints; ++i, ++position) {
      int segment{ bumps.getNearestSegmentIndex(*position) };
      auto bump{ bumps.getAtSegmentIndex(segment) };

      if (bump.inBump(*position)) {
        int currentCount{ bumpAccumulants.contains(segment) ? bumpAccumulants.at(segment) : 0 };
        bumpAccumulants.insert_or_assign(segment, currentCount + 1);
        ++totalInBumps;
      }
    }

    double entropy{ 0 };
    for (const auto& [_, numPointsInBump] : bumpAccumulants) {
      double probability{ static_cast<double>(numPointsInBump) / totalInBumps };
      entropy -= probability * std::log(probability);
    }

    return entropy;
  }

};
