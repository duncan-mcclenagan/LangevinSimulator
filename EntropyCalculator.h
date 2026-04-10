#pragma once

#include <unordered_map>
#include "PhaseCloud.h"
#include "BumpDistribution.h"

class EntropyCalculator {
private:
  std::unordered_map<int, int> accumulants{};

public:
  double calculateEntropy(const PhaseCloud& phaseCloud, const BumpDistribution& bumps, int particleIndex, int spatialCoordinateIndex);

};
