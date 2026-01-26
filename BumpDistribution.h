#pragma once

#include <array>
#include <algorithm>
#include <cmath>
#include <mdspan>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include "ArraySpan.h"
#include "Bump.h"
#include "ext/rapidcsv/src/rapidcsv.h"

template <int NSpatialCoords = 2>
requires (NSpatialCoords > 0)
class BumpDistribution {
private:
  // @todo Add inner const qualifiers (??)
  const std::vector<double> squareWidths;
  const std::vector<double> heights;
  const std::vector<std::array<double, NSpatialCoords>> positions;
  const std::unordered_map<ArraySpan<int, NSpatialCoords>, int> segmentToBumpId;

  const double inverseSegmentSpacing;
  const int numBumps;
  const int numParticles;
  const int numSpatialCoords;

  BumpDistribution(std::vector<double> squareWidths,
    std::vector<double> heights,
    std::vector<std::array<double, NSpatialCoords>> positions,
    double inverseSegmentSpacing,
    int numBumps,
    int numParticles,
    int numSpatialCoords)
    : squareWidths(squareWidths),
      heights(heights),
      positions(positions),
      inverseSegmentSpacing(inverseSegmentSpacing),
      numBumps(numBumps),
      numParticles(numParticles),
      numSpatialCoords(numSpatialCoords)
  {}

public:

  std::optional<int> getCandidateBumpId(const ArraySpan<double, NSpatialCoords> particlePositions) const {
    assert(particlePositions.size() == numParticles);

    static thread_local std::vector<std::array<int, NSpatialCoords>> segmentId{};
    if (segmentId.size() < particlePositions.size()) {
      segmentId.reserve(particlePositions.size());
    }

    for (int i{ 0 }; i < numParticles; ++i) {
      segmentId.push_back(std::array<int, NSpatialCoords>{});
      auto& discreteCoords{ segmentId[i] };
      auto& coords{ particlePositions[i] };
      for (int j{ 0 }; j < NSpatialCoords; ++j) {
        discreteCoords[j] = static_cast<int>(std::lround(coords[j] * inverseSegmentSpacing));
      }
    }

    ArraySpan<int, NSpatialCoords> segmentIdSpan{ segmentId.data(), numParticles };
    return segmentToBumpId.contains(segmentIdSpan)) ? segmentToBumpId.at(segmentIdSpan) : std::nullopt;
  }

  std::optional<int> getBumpId(const ArraySpan<double, NSpatialCoords> particlePositions) const {
    const auto candidateBumpId{ getBumpId(particlePositions) };
    if (candidateBumpId == std::nullopt) {
      return std::nullopt;
    }

    const double bumpSquareWidth{ squareWidths[candidateBumpId] };
    const ArraySpan<double, NSpatialCoords> bumpPosition{ &positions[candidateBumpId], numParticles };

    double squareDistance{ 0 };
    for (int i{ 0 }; i < numParticles; ++i) {
      for (int j{ 0 }; j < NSpatialCoords; ++j) {
        const double diff{ particlePositions[i, j] - bumpPosition[i, j] };
        squareDistance += diff * diff;
      }
    }

    return (squareDistance <= bumpSquareWidth) ? candidateBumpId : std::nullopt;
  }

  void computeForce(const ArraySpan<double, NSpatialCoords> particlePositions, ArraySpan<double, NSpatialCoords> forceOutput) const {
    assert(particlePositions.size() == forceOutput.size() == numParticles);

    // Start by getting the candidate bump.
    auto candidateBumpId{ getCandidateBumpId(particlePositions) };
    if (!candidateBumpId.has_value()) {
      for (int i{ 0 }; i < numParticles; ++i) {
        for (int j{ 0 }; j < NSpatialCoords; ++j) {
          forceOutput[i, j] = 0;
        }
      }
      return;
    }

    // Calculate the difference part of the force output first, and simultaneously
    // calculate the square distance. If we're outside the bump, return a zero force
    // vector.
    const double bumpSquareWidth{ squareWidths[candidateBumpId] };
    const ArraySpan<double, NSpatialCoords> bumpPosition{ &positions[candidateBumpId], numParticles };
    double squareDistance{ 0 };
    for (int i{ 0 }; i < numParticles; ++i) {
      for (int j{ 0 }; j < NSpatialCoords; ++j) {
        const double diff{ particlePositions[i, j] - bumpPosition[i, j] };
        squareDistance += diff * diff;
        forceOutput[i, j] = diff;
      }
    }
    if (squareDistance > bumpSquareWidth) {
      for (int i{ 0 }; i < numParticles; ++i) {
        for (int j{ 0 }; j < NSpatialCoords; ++j) {
          forceOutput[i, j] = 0;
        }
      }
      return;
    }

    // Use the bump function V(r) = -H*e^(1/((|r - a|/w)^2 - 1)) to produce our minimum.
    // `H` is the height of the bump, `w` is the width, and `a` is the center.
    // The corresponding force vector is -H*e^(1/(|r - a|^2/w^2 - 1))*(2*(x - a_x)/w^2)/(|r - a|^2/w^2 - 1)^2
    // for each component `x`.
    const double inverseBumpSquareWidth{ 1 / bumpSquareWidth };
    const double expArg{ 1 / (squareDistance * inverseBumpSquareWidth - 1) };
    const double bumpHeight{ heights[candidateBumpId] };

    const double multiplier{ -bumpHeight * std::exp(expArg) * (2 * inverseBumpSquareWidth) * expArg * expArg };
    for (int i{ 0 }; i < numParticles; ++i) {
      for (int j{ 0 }; j < NSpatialCoords; ++j) {
        forceOutput[i, j] *= multiplier;
      }
    }
  }

  static BumpDistribution fromFile(const char* filename, double segmentSpacing) {
    rapidcsv::Document csv(filename, rapidcsv::LabelParams(), rapidcsv::SeparatorParams(), rapidcsv::LineReaderParams(pSkipEmptyLines = true));

    int widthIndex{ csv.GetColumnIdx("width") };
    int centerIndex{ csv.GetColumnIdx("center") };
    if (widthIndex < 0 || centerIndex < 0) {
      throw std::runtime_error("CSV file must contain 'width' and 'center' columns.");
    }

    std::vector<int> positionIndices{};

    // Determine the dimension by reading the column headers.
    int currentSpatialCoordinate{ 0 };
    int currentParticleIndex{ 0 };
    std::optional<int> spatialDimension{};
    int columnIndex{ csv.GetColumnIdx("pos-0-0") };
    if (columnIndex < 0) {
      throw std::runtime_error("CSV file must contain at least the first position column (pos-0-0).");
    }
    do {
      // @todo Something a bit cleaner here?
      positionIndices.push_back(columnIndex);
      columnIndex = csv.GetColumnIdx(getCsvCoordinateColumnName(currentParticleIndex, ++currentSpatialCoordinate));
      if (columnIndex < 0) {
        if (spatialDimension.has_value()) {
          if (currentSpatialCoordinate != spatialDimension) {
            throw std::runtime_error("Inconsistent number of spatial dimensions in position columns.");
          }
        }
        else spatialDimension = currentSpatialCoordinate;

        currentSpatialCoordinate = 0;
        columnIndex = csv.GetColumnIdx(getCsvCoordinateColumnName(++currentParticleIndex, currentSpatialCoordinate));
      }
    } while (columnIndex >= 0);
    int numParticles{ currentParticleIndex };
    assert(spatialDimension.has_value());

    int numBumps{ csv.GetRowCount() };
    // @todo Proper trapping here...
    std::vector<double> widths{ csv.GetColumn<double>("width") };
    std::vector<double> heights{ csv.GetColumn<double>("height") };
    if (!(widths.size() == heights.size() == numBumps)) {
      throw std::runtime_error("Inconsistent or invalid number of entries in width/height columns.");
    }

    // @todo Something cleaner here??
    std::vector<std::array<double, NSpatialCoords>> positions{};
    for (int i{ 0 }; i < numBumps; ++i) {
      for (int j{ 0 }; j < numParticles; ++j) {
        for (int k{ 0 }; k < spatialDimension; ++k)
      }
    }
  }

private:
  static std::string getCsvCoordinateColumnName(int particleIndex, int spatialCoordinateIndex) {
    return "pos-" + std::to_string(particleIndex) + "-" + std::to_string(spatialCoordinateIndex);
  }

};
