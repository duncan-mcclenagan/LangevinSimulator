#include <assert.h>
#include <optional>
#include <stdexcept>
#include "ext/rapidcsv/src/rapidcsv.h"
#include "BumpDistribution.h"

void BumpDistribution::computeForce(CoordinateSpan<const double> particlePositions, CoordinateSpan<double> forceOutput) const {
  assert(particlePositions.totalSize() == forceOutput.totalSize() == totalPhaseSpaceDimension);

  double* forceOutputData{ forceOutput.raw() };
  const double* particlePositionsData{ particlePositions.raw() };

  auto bumpResearchResult{ bumpsByPosition.find(particlePositions.raw()) };

  if (bumpResearchResult == bumpsByPosition.end()) {
    // Nothing was found => no force.
    for (int i{ 0 }; i < totalPhaseSpaceDimension; ++i) {
      forceOutputData[i] = 0;
    }
    return;
  }

  const double* bumpPosition{ bumpResearchResult->first };
  const auto bumpData{ bumpResearchResult->second };

  // Calculate the difference part of the force output first, and simultaneously
  // calculate the square distance. If we're outside the bump, return a zero force
  // vector.
  const double bumpSquareWidth{ bumpData.squareWidth };
  double squareDistance{ 0 };
  for (int i{ 0 }; i < totalPhaseSpaceDimension; ++i) {
    const double diff{ particlePositionsData[i] - bumpPosition[i] };
    squareDistance += diff * diff;
    forceOutputData[i] = diff;
  }

  if (squareDistance > bumpSquareWidth) {
    for (int i{ 0 }; i < totalPhaseSpaceDimension; ++i) {
      forceOutputData[i] = 0;
    }
    return;
  }

  // Use the bump function V(r) = -H*e^(1/((|r - a|/w)^2 - 1)) to produce our minimum.
  // `H` is the height of the bump, `w` is the width, and `a` is the center.
  // The corresponding force vector is -H*e^(1/(|r - a|^2/w^2 - 1))*(2*(x - a_x)/w^2)/(|r - a|^2/w^2 - 1)^2
  // for each component `x`.
  const double inverseBumpSquareWidth{ 1 / bumpSquareWidth };
  const double expArg{ 1 / (squareDistance * inverseBumpSquareWidth - 1) };
  const double bumpHeight{ bumpData.height };

  const double multiplier{ -bumpHeight * std::exp(expArg) * (2 * inverseBumpSquareWidth) * expArg * expArg };
  for (int i{ 0 }; i < totalPhaseSpaceDimension; ++i) {
    forceOutputData[i] *= multiplier;
  }
}

std::generator<bool> BumpDistribution::enumerateDiscreteCoordinates(CoordinateSpan<int> discreteCoordinates) const {
  assert(discreteCoordinates.totalSize() == totalPhaseSpaceDimension);

  for (int i{ 0 }; i < numBumps; ++i) {
    computeDiscreteCoordinates(CoordinateSpan<const double>(positions.data() + i * totalPhaseSpaceDimension, numParticles, numSpatialCoords), discreteCoordinates);
    co_yield true;
  }
}

void BumpDistribution::computeDiscreteCoordinates(CoordinateSpan<const double> coordinates, CoordinateSpan<int> discreteCoordinates) const {
  const int totalSize{ coordinates.totalSize() };
  assert(totalSize == discreteCoordinates.totalSize());

  const double* coordinatesData{ coordinates.raw() };
  int* discreteCoordinatesData{ discreteCoordinates.raw() };
  for (int i{ 0 }; i < totalSize; ++i) {
    discreteCoordinatesData[i] = toDiscrete(coordinatesData[i], inverseSegmentSpacing);
  }
}

void BumpDistribution::toFile(const char* filename) const {
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
  for (const auto& [rawPosition, bumpData] : bumpsByPosition) {
    csv.SetCell(widthIndex, i, std::sqrt(bumpData.squareWidth));
    csv.SetCell(heightIndex, i, bumpData.height);
    CoordinateSpan<const double> position{ rawPosition, numParticles, numSpatialCoords };
    for (int j{ 0 }, n{ startPositionIndex }; j < numParticles; ++j) {
      for (int k{ 0 }; k < numSpatialCoords; ++k, ++n) {
        csv.SetCell(n, i, position[j, k]);
      }
    }

    ++i;
  }
}

BumpDistribution BumpDistribution::fromFile(const char* filename, double segmentSpacing) {
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
  const double inverseSegmentSpacing{ 1 / segmentSpacing };
  const int totalPhaseSpaceDimension{ spatialDimension.value() * numParticles };
  BumpMap bumpsByPosition(numBumps, BumpDataHash(totalPhaseSpaceDimension, inverseSegmentSpacing), BumpDataEquality(totalPhaseSpaceDimension, inverseSegmentSpacing));
  bumpsByPosition.reserve(numBumps);
  positions.reserve(totalPhaseSpaceDimension * numBumps);
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

    if (bumpsByPosition.contains(positionStart)) {
      throw std::runtime_error("Input CSV contains multiple bumps in one segment.");
    }
    bumpsByPosition.insert({ positionStart, BumpData(width, height) });
    return BumpDistribution(std::move(positions), std::move(bumpsByPosition), inverseSegmentSpacing, numBumps, numParticles, spatialDimension.value());
  }
}

inline bool BumpDistribution::BumpDataEquality::operator()(const double* firstPosition, const double* secondPosition) const {
  for (int i{ 0 }; i < totalPhaseSpaceDimension; ++i) {
    if (BumpDistribution::toDiscrete(firstPosition[i], inverseSegmentSpacing) != BumpDistribution::toDiscrete(secondPosition[i], inverseSegmentSpacing)) {
      return false;
    }
  }

  return true;
}

inline size_t BumpDistribution::BumpDataHash::operator()(const double* position) const {
  size_t hash{ 0 };
  for (int i{ 0 }; i < totalPhaseSpaceDimension; ++i) {
    hash ^= std::hash<int>{}(BumpDistribution::toDiscrete(position[i], inverseSegmentSpacing)) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
  }
  return hash;
}
