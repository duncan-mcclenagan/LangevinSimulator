#include <assert.h>
#include <cmath>
#include "EntropyCalculator.h"

double EntropyCalculator::calculateEntropy(const PhaseCloud& phaseCloud, const BumpDistribution& bumps, int particleIndex, int spatialCoordinateIndex) {
  assert(phaseCloud.numberOfParticles == bumps.numParticles);
  assert(phaseCloud.spatialDimension == bumps.numSpatialCoords);
  assert(particleIndex >= 0 && particleIndex < phaseCloud.numberOfParticles);
  assert(spatialCoordinateIndex >= 0 && spatialCoordinateIndex < phaseCloud.spatialDimension);

  accumulants.clear();

  for (auto position : phaseCloud.enumeratePointPositionsConst()) {
    int coordinate{ bumps.toDiscrete(position[particleIndex, spatialCoordinateIndex]) };
    int currentCount{ accumulants.contains(coordinate) ? accumulants.at(coordinate) : 0 };
    accumulants.insert_or_assign(coordinate, currentCount + 1);
  }

  double entropy{ 0 };
  double totalNumPoints{ static_cast<double>(phaseCloud.numberOfPoints) };
  for (const auto& [_, numPoints] : accumulants) {
    double probability{ numPoints / totalNumPoints };
    entropy -= probability * std::log(probability);
  }

  return entropy;
}
