#pragma once

#include <array>
#include <algorithm>
#include <cmath>
#include <mdspan>
#include <optional>
#include <stdexcept>
#include <generator>
#include <string>
#include <unordered_set>
#include "CoordinateSpan.h"
#include "Bump.h"
#include "ext/rapidcsv/src/rapidcsv.h"

class BumpDistribution {
private:
  struct BumpData {
    const double squareWidth;
    const double height;

    constexpr BumpData(double width, double height)
      : squareWidth(width * width), height(height > 0 ? height : -height) {}

  };

  class BumpDataEquality {

  private:
    const int totalPhaseSpaceDimension;
    const double inverseSegmentSpacing;

  public:
    BumpDataEquality(int totalPhaseSpaceDimension, double inverseSegmentSpacing)
      : totalPhaseSpaceDimension(totalPhaseSpaceDimension),
        inverseSegmentSpacing(inverseSegmentSpacing)
    {}

    bool operator()(const double* firstPosition, const double* secondPosition) const;

  };

  class BumpDataHash {

  private:
    const int totalPhaseSpaceDimension;
    const double inverseSegmentSpacing;

  public:
    BumpDataHash(int totalPhaseSpaceDimension, double inverseSegmentSpacing)
      : totalPhaseSpaceDimension(totalPhaseSpaceDimension),
        inverseSegmentSpacing(inverseSegmentSpacing)
    {}

    size_t operator()(const double* position) const;

  };

public:
  const int numBumps;
  const int numParticles;
  const int numSpatialCoords;
  const int totalPhaseSpaceDimension;
private:
  const std::vector<double> positions;

  typedef std::unordered_map<const double*, BumpData, BumpDataHash, BumpDataEquality> BumpMap;
  const BumpMap bumpsByPosition;

  const double inverseSegmentSpacing;

  BumpDistribution(std::vector<double> positions,
    BumpMap bumpsByPosition,
    double inverseSegmentSpacing,
    int numBumps,
    int numParticles,
    int numSpatialCoords)
    : positions(positions),
      bumpsByPosition(bumpsByPosition),
      inverseSegmentSpacing(inverseSegmentSpacing),
      numBumps(numBumps),
      numParticles(numParticles),
      numSpatialCoords(numSpatialCoords),
      totalPhaseSpaceDimension(numParticles * numSpatialCoords)
  {}

public:

  void computeForce(CoordinateSpan<const double> particlePositions, CoordinateSpan<double> forceOutput) const;

  std::generator<bool> enumerateDiscreteCoordinates(CoordinateSpan<int> discreteCoordinates) const;

  void toFile(const char* filename) const;

  template<typename F>
    requires std::invocable<F, double&, double&, CoordinateSpan<double>>
  static BumpDistribution generate(int numParticles, int spatialDimension, double segmentSpacing, int numSegmentsInOneDimension, double bumpProbability, F bumpGenerator) {
    // @todo Argument error handling.
    const int totalPhaseSpaceDimension{ numParticles * spatialDimension };
    const int estimatedNumBumps{ static_cast<int>(std::pow(numSegmentsInOneDimension * bumpProbability, totalPhaseSpaceDimension)) };
    const double inverseSegmentSpacing{ 1 / segmentSpacing };

    std::vector<double> positions{};
    BumpMap bumpsByPosition(totalPhaseSpaceDimension, BumpDataHash(totalPhaseSpaceDimension, inverseSegmentSpacing), BumpDataEquality(totalPhaseSpaceDimension, inverseSegmentSpacing));
    positions.reserve(estimatedNumBumps * totalPhaseSpaceDimension);

    std::mt19937_64 randomGenerator((std::random_device())());
    std::bernoulli_distribution bumpDecider(bumpProbability);

    // @todo Keep track of "current base position" and run through the entire grid, updating the rightmost base position
    // value that does not push us beyond the range (otherwise reset it).
  }

  static BumpDistribution fromFile(const char* filename, double segmentSpacing);

private:
  void computeDiscreteCoordinates(CoordinateSpan<const double> coordinates, CoordinateSpan<int> discreteCoordinates) const;

  static int toDiscrete(double value, double inverseSegmentSpacing) {
    return static_cast<int>(std::floor(value * inverseSegmentSpacing));
  }

  static std::string getCsvCoordinateColumnName(int particleIndex, int spatialCoordinateIndex) {
    return "pos-" + std::to_string(particleIndex) + "-" + std::to_string(spatialCoordinateIndex);
  }

};
