#pragma once

#include <array>
#include <algorithm>
#include <cmath>
#include <mdspan>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include "CoordinateSpan.h"
#include "Bump.h"
#include "ext/rapidcsv/src/rapidcsv.h"

class BumpDistribution {
public:
  const int numBumps;
  const int numParticles;
  const int numSpatialCoords;
private:
  const std::vector<double> positions;
  const std::unordered_map<const double*, BumpData> bumpsByPosition;

  const double inverseSegmentSpacing;

  BumpDistribution(std::vector<double> positions,
    std::unordered_map<const double*, BumpData> bumpsByPosition,
    double inverseSegmentSpacing,
    int numBumps,
    int numParticles,
    int numSpatialCoords)
    : squareWidths(squareWidths),
      heights(heights),
      positions(positions),
      discretePositions(discretePositions),
      bumpsByPosition(bumpsByPosition),
      inverseSegmentSpacing(inverseSegmentSpacing),
      numBumps(numBumps),
      numParticles(numParticles),
      numSpatialCoords(numSpatialCoords)
  {}

public:

  std::optional<BumpData&> getBumpInSegment(const CoordinateSpan<double> particlePositions) const {
    assert(particlePositions.size() == numParticles);

    static thread_local std::vector<int> segmentId{};
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

    CoordinateSpan<int, NSpatialCoords> segmentIdSpan{ segmentId.data(), numParticles };
    return bumpsByPosition.contains(segmentIdSpan)) ? bumpsByPosition.at(segmentIdSpan) : std::nullopt;
  }

  std::optional<BumpData&> getBump(const CoordinateSpan<double, NSpatialCoords> particlePositions) const {
    const auto candidateBumpId{ getBumpId(particlePositions) };
    if (candidateBumpId == std::nullopt) {
      return std::nullopt;
    }

    const double bumpSquareWidth{ squareWidths[candidateBumpId] };
    const CoordinateSpan<double, NSpatialCoords> bumpPosition{ &positions[candidateBumpId], numParticles };

    double squareDistance{ 0 };
    for (int i{ 0 }; i < numParticles; ++i) {
      for (int j{ 0 }; j < NSpatialCoords; ++j) {
        const double diff{ particlePositions[i, j] - bumpPosition[i, j] };
        squareDistance += diff * diff;
      }
    }

    return (squareDistance <= bumpSquareWidth) ? candidateBumpId : std::nullopt;
  }

  void computeForce(const CoordinateSpan<double, NSpatialCoords> particlePositions, CoordinateSpan<double, NSpatialCoords> forceOutput) const {
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
    const CoordinateSpan<double, NSpatialCoords> bumpPosition{ &positions[candidateBumpId], numParticles };
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

  void toFile(const char* filename) const {
    // @todo Proper error handling...
    rapidcsv::Document csv(filename);
    csv.Clear();

    const int widthIndex{ 0 };
    const int heightIndex{ 1 };
    const int startPositionIndex{ 2 };
    csv.SetColumnName(widthIndex, "width");
    csv.SetColumnName(heightIndex, "height");
    for (int i{ 0 }, n{ startPositionIndex }; i < numParticles; ++i) {
       for (int j{ 0 }; j < numSpatialCoords; ++j, ++n) {
         csv.SetColumnName(n, getCsvCoordinateColumnName(i, j));
       }
    }

    int i{ 0 };
    for (const auto& [_, bumpData] : bumpsByPosition) {
      csv.SetCell(widthIndex, i, std::sqrt(bumpData.squareWidth));
      csv.SetCell(heightIndex, i, bumpData.height);
      CoordinateSpan<const double> position{ bumpData.center, numParticles, numSpatialCoords };
      for (int j{ 0 }, n{ startPositionIndex }; j < numParticles; ++j) {
        for (int k{ 0 }; k < numSpatialCoords; ++k, ++n) {
          csv.SetCell(n, i, position[j, k]);
        }
      }

      ++i;
    }
  }

  template<std::invocable<bool, double&, double&, CoordinateSpan<double>> BumpGenerator>
  static BumpDistribution generate(int numParticles, int spatialDimension, double segmentSpacing, BumpGenerator generator) {
    constexpr int defaultApprxNumBumps{ 100 };
    // @todo Argument error handling.
    std::vector<double> positions{};
    std::unordered_map<const double*, BumpData> bumpsByPosition{};
    const int coordinatesPerBump{ numParticles * spatialDimension };
    positions.reserve(coordinatesPerBump * defaultApprxNumBumps);
    positions.resize(coordinatesPerBump, 1);

    double width{ 1 }, height{ 1 };
    
    int i{ 0 };
    bool generateNewBump{ true };
    while (generateNewBump) {
      positions.resize(coordinatesPerBump);
      double* positionStart{ &positions[i * coordinatesPerBump] };
      CoordinateSpan<double> position(positionStart, numParticles, spatialDimension);
      generateNewBump = generator(width, height, position);
      if (bumpsByPosition.contains(positionStart)) {
        throw std::runtime_error("Only one bump per segment permitted.");
      }

      bumpsByPosition.insert({ positionStart, BumpData(positionStart, width, height) });

      ++i;
    }

    return BumpDistribution();
  } 

  static BumpDistribution fromFile(const char* filename, double segmentSpacing) {
    rapidcsv::Document csv(filename,
      rapidcsv::LabelParams(),
      rapidcsv::SeparatorParams(),
      rapidcsv::ConverterParams(),
      rapidcsv::LineReaderParams(false, false, true));

    int widthIndex{ csv.GetColumnIdx("width") };
    int heightIndex{ csv.GetColumnIdx("height") };
    if (widthIndex < 0 || heightIndex < 0) {
      throw std::runtime_error("CSV file must contain 'width' and 'height' columns.");
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
    // @todo Something cleaner here??
    // @todo Proper trapping with invalid conversions during GetCell().
    std::vector<double> positions{};
    std::unordered_map<const double*, BumpData> bumpsByPosition{};
    bumpsByPosition.reserve(numBumps);
    positions.reserve(spatialDimension.value() * numParticles * numBumps);
    for (int i{ 0 }; i < numBumps; ++i) {
      double width{ csv.GetCell<double>(widthIndex, i) };
      double height{ csv.GetCell<double>(height, i)};

      double* positionStart{ &positions[positions.size()] };
      for (int j{ 0 }, n{ 0 }; j < numParticles; ++j) {
        for (int k{ 0 }; k < spatialDimension; ++k, ++n) {
          double positionComponent{ csv.GetCell<double>(positionIndices[n], i) };
          positions.push_back(positionComponent);
        }
      }

      CoordinateSpan<const double> position(positionStart, numParticles, spatialDimension.value());
      if (bumpsByPosition.contains(position)) {
        throw std::runtime_error("Input CSV contains multiple bumps in one segment.");
      }
      bumpsByPosition.insert({ position, Bump(position, width, height) });

      return BumpDistribution();
    }


  }

private:
  void computeDiscreteCoordinates(CoordinateSpan<const double> coordinates, CoordinateSpan<int> discreteCoordinates) {

  }

  static std::string getCsvCoordinateColumnName(int particleIndex, int spatialCoordinateIndex) {
    return "pos-" + std::to_string(particleIndex) + "-" + std::to_string(spatialCoordinateIndex);
  }

};
