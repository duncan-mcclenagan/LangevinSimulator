#pragma once

#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include "PhaseCloud.h"
#include "BumpDistribution.h"

class EntropyCalculator {
private:
  std::unordered_map<void*, int> bumpAccumulants{};

public:
  double calculateEntropy(const PhaseCloud& phaseCloud, const BumpDistribution& bumps) {
    assert(phaseCloud.numberOfParticles == bumps.numParticles);
    assert(phaseCloud.spatialDimension == bumps.numSpatialCoords);

    bumpAccumulants.clear();

    int totalInBumps{ 0 };
    for (auto position : phaseCloud.enumeratePointPositionsConst()) {
      void* bump{ bumps.getBumpIdentifier(position) };
      if (bump != nullptr) {
        int currentCount{ bumpAccumulants.contains(bump) ? bumpAccumulants.at(bump) : 0 };
        bumpAccumulants.insert_or_assign(bump, currentCount + 1);
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
