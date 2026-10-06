#pragma once

#include <generator>
#include <string>
#include <random>
#include <unordered_set>
#include <vector>
#include <concepts>
#include <algorithm>
#include <cmath>
#include <utility>
#include "CoordinateSpan.h"

// @todo-H The bumpsByPosition thing doesn't work as-is because it doesn't reference stable memory locations.
// Instead, we should simply use bump indices and access the bump data from the positions vector.

class BumpDistribution {
private:
  struct BumpData {
    const double squareWidth;
    const double height;
    const int startPositionIndex;

    constexpr BumpData(double width, double height, int startPositionIndex)
      : squareWidth(width * width), height(height > 0 ? height : -height), startPositionIndex(startPositionIndex >= 0 ? startPositionIndex : 0) {}

  };

  class BumpDataEquality {

    // Enable heterogenous lookup.
    using is_transparent = void;

  private:
    bool inSameSector(CoordinateSpan<const double> position, const BumpData& data, const BumpDistribution& distribution) const;

  public:

    bool operator()(const BumpData& first, const BumpData& second) const;

    inline bool operator()(const BumpData& data, std::pair<CoordinateSpan<const double>, const BumpDistribution&> positionAndDistribution) const {
      return inSameSector(positionAndDistribution.first, data, positionAndDistribution.second);
    }

    inline bool operator()(std::pair<CoordinateSpan<const double>, const BumpDistribution&> positionAndDistribution, const BumpData& data) const {
      return inSameSector(positionAndDistribution.first, data, positionAndDistribution.second);
    }

    inline bool operator()(const BumpData& data, std::pair<CoordinateSpan<double>, const BumpDistribution&> positionAndDistribution) const {
      return inSameSector(positionAndDistribution.first.asConst(), data, positionAndDistribution.second);
    }

    inline bool operator()(std::pair<CoordinateSpan<double>, const BumpDistribution&> positionAndDistribution, const BumpData& data) const {
      return inSameSector(positionAndDistribution.first.asConst(), data, positionAndDistribution.second);
    }

  };

  class BumpDataHash {
    using is_transparent = void;

  private:
    size_t getHashForPosition(CoordinateSpan<const double> position, double inverseSegmentSpacing, int totalPhaseSpaceDimension) const;
    inline size_t getHashForPosition(CoordinateSpan<const double> position, const BumpDistribution& distribution) const {
      return getHashForPosition(position, distribution.inverseSegmentSpacing, distribution.totalPhaseSpaceDimension);
    }

  public:
    inline size_t operator()(const BumpData& data) const {
      data.
    }

    inline size_t operator()(std::pair<CoordinateSpan<const double>, const BumpDistribution&> positionAndDistribution) const {
      return getHashForPosition(positionAndDistribution.first, positionAndDistribution.second);
    }

    size_t operator()(std::pair<CoordinateSpan<double>, const BumpDistribution&> positionAndDistribution) const {
      return getHashForPosition(positionAndDistribution.first.asConst(), positionAndDistribution.second);
    }

  };

public:
  const int numBumps;
  const int numParticles;
  const int numSpatialCoords;
  const int totalPhaseSpaceDimension;
private:
  const std::vector<double> positions;

  typedef std::unordered_set<BumpData, BumpDataHash, BumpDataEquality> BumpSet;
  const BumpSet bumps;

  const double inverseSegmentSpacing;

  BumpDistribution(std::vector<double> positions,
    BumpSet bumps,
    double inverseSegmentSpacing,
    int numBumps,
    int numParticles,
    int numSpatialCoords)
    : positions(positions),
      bumps(bumps),
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
    // Of form (double &width, double &height, CoordinateSpan<double> position)
    requires std::invocable<F, double&, double&, CoordinateSpan<double>>
  static BumpDistribution generate(int numParticles, int spatialDimension, double segmentSpacing, int numSegmentsInOneDimension, double bumpProbability, F bumpGenerator) {
    // @todo Argument error handling.
    const int totalPhaseSpaceDimension{ numParticles * spatialDimension };
    const int estimatedNumBumps{ static_cast<int>(std::pow(numSegmentsInOneDimension * bumpProbability, totalPhaseSpaceDimension)) };
    const double inverseSegmentSpacing{ 1 / segmentSpacing };

    std::vector<double> positions{};
    BumpSet bumps(estimatedNumBumps);
    positions.reserve(static_cast<size_t>(estimatedNumBumps) * totalPhaseSpaceDimension);

    std::mt19937_64 randomGenerator((std::random_device())());
    std::bernoulli_distribution bumpDecider(bumpProbability);

    std::vector<int> indices(totalPhaseSpaceDimension);
    while (indices[0] < numSegmentsInOneDimension) {
      if (bumpDecider(randomGenerator)) {
        positions.reserve(positions.size() + totalPhaseSpaceDimension);
        for (int i{ 0 }; i < totalPhaseSpaceDimension; ++i) {
          positions.push_back(indices[i] * segmentSpacing);
        }

        double bumpWidth{};
        double bumpHeight{};
        int bumpPositionStartIndex{ static_cast<int>(positions.size()) - totalPhaseSpaceDimension };
        CoordinateSpan<double> bumpPosition(positions.data() + bumpPositionStartIndex, numParticles, spatialDimension);
        bumpGenerator(bumpWidth, bumpHeight, bumpPosition);

        bumps.insert(BumpData(bumpWidth, bumpHeight, bumpPositionStartIndex));
      }

      // Advance right-most position unless it takes us beyond the range, in which case reset it and advance the next right-most position, etc.
      for (int i{ totalPhaseSpaceDimension - 1 }; i >= 0; --i) {
        if (++indices[i] < numSegmentsInOneDimension) break;
        else indices[i] = 0;
      }
    }

    return BumpDistribution(std::move(positions),
      bumps,
      inverseSegmentSpacing,
      static_cast<int>(bumps),
      numParticles,
      spatialDimension
    );
  }

  static BumpDistribution fromFile(const char* filename, double segmentSpacing);

  int toDiscrete(double value) const {
    return toDiscrete(value, inverseSegmentSpacing);
  }

private:
  void computeDiscreteCoordinates(CoordinateSpan<const double> coordinates, CoordinateSpan<int> discreteCoordinates) const;

  static int toDiscrete(double value, double inverseSegmentSpacing) {
    return static_cast<int>(std::floor(value * inverseSegmentSpacing));
  }

  static std::string getCsvCoordinateColumnName(int particleIndex, int spatialCoordinateIndex) {
    return "pos-" + std::to_string(particleIndex) + "-" + std::to_string(spatialCoordinateIndex);
  }

};
