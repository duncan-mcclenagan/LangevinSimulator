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
#include <memory>
#include "CoordinateSpan.h"

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
  private:
    const std::shared_ptr<const std::vector<double>> positions;
    const double inverseSegmentSpacing;
    const int totalPhaseSpaceDimension;

    bool inSameSector(CoordinateSpan<const double> position, const BumpData& data) const;

  public:
    using is_transparent = void;

    BumpDataEquality(std::shared_ptr<const std::vector<double>> positions, double inverseSegmentSpacing, int totalPhaseSpaceDimension)
      : positions(positions), inverseSegmentSpacing(inverseSegmentSpacing), totalPhaseSpaceDimension(totalPhaseSpaceDimension) {}

    bool operator()(const BumpData& first, const BumpData& second) const;

    bool operator()(const BumpData& data, CoordinateSpan<const double> position) const {
      return inSameSector(position, data);
    }

    bool operator()(CoordinateSpan<const double> position, const BumpData& data) const {
      return inSameSector(position, data);
    }

    bool operator()(const BumpData& data, CoordinateSpan<double> position) const {
      return inSameSector(position.asConst(), data);
    }

    bool operator()(CoordinateSpan<double> position, const BumpData& data) const {
      return inSameSector(position.asConst(), data);
    }

  };

  class BumpDataHasher {
  private:
    const std::shared_ptr<const std::vector<double>> positions;
    const double inverseSegmentSpacing;
    const int totalPhaseSpaceDimension;

    size_t getHashForPosition(CoordinateSpan<const double> position) const;

  public:
    using is_transparent = void;

    BumpDataHasher(std::shared_ptr<const std::vector<double>> positions, double inverseSegmentSpacing, int totalPhaseSpaceDimension)
      : positions(positions), inverseSegmentSpacing(inverseSegmentSpacing), totalPhaseSpaceDimension(totalPhaseSpaceDimension) {}

    size_t operator()(const BumpData& data) const;

    size_t operator()(CoordinateSpan<const double> position) const {
      return getHashForPosition(position);
    }

    size_t operator()(CoordinateSpan<double> position) const {
      return getHashForPosition(position.asConst());
    }

  };

private:
  const std::shared_ptr<const std::vector<double>> positions;
  const double inverseSegmentSpacing;

  typedef std::unordered_set<BumpData, BumpDataHasher, BumpDataEquality> BumpSet;
  const BumpSet bumps;

public:
  const int numBumps;
  const int numParticles;
  const int numSpatialCoords;
  const int totalPhaseSpaceDimension;

private:
  BumpDistribution(std::shared_ptr<const std::vector<double>> positions,
    BumpSet bumps,
    double inverseSegmentSpacing,
    int numParticles,
    int numSpatialCoords)
    : positions(std::move(positions)),
      inverseSegmentSpacing(inverseSegmentSpacing),
      bumps(std::move(bumps)),
      numBumps(static_cast<int>(this->bumps.size())),
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
    const int estimatedNumBumps{ static_cast<int>(std::pow(numSegmentsInOneDimension, totalPhaseSpaceDimension) * bumpProbability) };
    const double inverseSegmentSpacing{ 1 / segmentSpacing };

    auto positions{ std::make_shared<std::vector<double>>() };
    BumpSet bumps(estimatedNumBumps,
      BumpDataHasher(positions, inverseSegmentSpacing, totalPhaseSpaceDimension),
      BumpDataEquality(positions, inverseSegmentSpacing, totalPhaseSpaceDimension));
    positions->reserve(static_cast<size_t>(estimatedNumBumps) * totalPhaseSpaceDimension);

    std::mt19937_64 randomGenerator((std::random_device())());
    std::bernoulli_distribution bumpDecider(bumpProbability);

    std::vector<int> indices(totalPhaseSpaceDimension);
    while (indices[0] < numSegmentsInOneDimension) {
      if (bumpDecider(randomGenerator)) {
        positions->reserve(positions->size() + totalPhaseSpaceDimension);
        for (int i{ 0 }; i < totalPhaseSpaceDimension; ++i) {
          positions->push_back(indices[i] * segmentSpacing);
        }

        double bumpWidth{};
        double bumpHeight{};
        int bumpPositionStartIndex{ static_cast<int>(positions->size()) - totalPhaseSpaceDimension };
        CoordinateSpan<double> bumpPosition(positions->data() + bumpPositionStartIndex, numParticles, spatialDimension);
        bumpGenerator(bumpWidth, bumpHeight, bumpPosition);

        bumps.insert(BumpData(bumpWidth, bumpHeight, bumpPositionStartIndex));
      }

      // Advance right-most position unless it takes us beyond the range, in which case reset it and advance the next right-most position, etc.
      // But don't reset the left index, because it needs to overflow to trigger the end of the loop.
      for (int i{ totalPhaseSpaceDimension - 1 }; i >= 0; --i) {
        if (++indices[i] < numSegmentsInOneDimension) break;
        else if (i > 0) indices[i] = 0;
      }
    }

    return BumpDistribution(std::move(positions),
      std::move(bumps),
      inverseSegmentSpacing,
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
