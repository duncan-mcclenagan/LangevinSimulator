#include <cassert>
#include <optional>
#include <stdexcept>
#include <cmath>
#include <generator>
#include <vector>
#include <utility>
#include <type_traits>
#include <functional>
#include "ext/rapidcsv/src/rapidcsv.h"
#include "CoordinateSpan.h"
#include "BumpDistribution.h"

void BumpDistribution::computeForce(CoordinateSpan<const double> particlePositions, CoordinateSpan<double> forceOutput) const {
  assert(particlePositions.totalSize() == forceOutput.totalSize() && forceOutput.totalSize() == totalPhaseSpaceDimension);

  double* forceOutputData{ forceOutput.raw() };
  const double* particlePositionsData{ particlePositions.raw() };

  auto bumpSearchResult{ bumps.find(particlePositions) };

  if (bumpSearchResult == bumps.end()) {
    // Nothing was found => no force.
    for (int i{ 0 }; i < totalPhaseSpaceDimension; ++i) {
      forceOutputData[i] = 0;
    }
    return;
  }

  const BumpData& bumpData{ *bumpSearchResult };
  const double* bumpPosition{ positions->data() + bumpData.startPositionIndex };

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

  for (const BumpData& bump : bumps) {
    computeDiscreteCoordinates(CoordinateSpan<const double>(positions->data() + bump.startPositionIndex, numParticles, numSpatialCoords), discreteCoordinates);
    co_yield true;
  }
}

void BumpDistribution::computeDiscreteCoordinates(CoordinateSpan<const double> coordinates, CoordinateSpan<int> discreteCoordinates) const {
  const int totalSize{ static_cast<int>(coordinates.totalSize()) };
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
  for (const auto& bumpData : bumps) {
    csv.SetCell(widthIndex, i, std::sqrt(bumpData.squareWidth));
    csv.SetCell(heightIndex, i, bumpData.height);
    CoordinateSpan<const double> position{ positions->data() + bumpData.startPositionIndex, numParticles, numSpatialCoords};
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

  int numBumps{ static_cast<int>(csv.GetRowCount()) };
  // @todo Something cleaner here??
  // @todo Proper trapping with invalid conversions during GetCell().
  auto positions{ std::make_shared<std::vector<double>>() };
  const double inverseSegmentSpacing{ 1 / segmentSpacing };
  const int totalPhaseSpaceDimension{ spatialDimension.value() * numParticles };
  BumpSet bumps(numBumps,
    BumpDataHasher(positions, inverseSegmentSpacing, totalPhaseSpaceDimension),
    BumpDataEquality(positions, inverseSegmentSpacing, totalPhaseSpaceDimension));
  positions->reserve(static_cast<size_t>(totalPhaseSpaceDimension) * numBumps);
  for (int i{ 0 }; i < numBumps; ++i) {
    double width{ csv.GetCell<double>(widthIndex, i) };
    double height{ csv.GetCell<double>(heightIndex, i)};

    int startPositionIndex{ static_cast<int>(positions->size()) };
    for (int j{ 0 }, n{ 0 }; j < numParticles; ++j) {
      for (int k{ 0 }; k < spatialDimension; ++k, ++n) {
        double positionComponent{ csv.GetCell<double>(positionIndices[n], i) };
        positions->push_back(positionComponent);
      }
    }

    double* positionStart{ positions->data() + startPositionIndex };
    CoordinateSpan<const double> position(positionStart, numParticles, spatialDimension.value());
    if (bumps.contains(position)) {
      throw std::runtime_error("Input CSV contains multiple bumps in one segment.");
    }
    bumps.insert(BumpData(width, height, startPositionIndex));
  }

  return BumpDistribution(positions, std::move(bumps), inverseSegmentSpacing, numParticles, spatialDimension.value());
}

bool BumpDistribution::BumpDataEquality::inSameSector(CoordinateSpan<const double> position, const BumpData& data) const {
  if (position.totalSize() != static_cast<size_t>(totalPhaseSpaceDimension)) {
    return false;
  }

  const double* const storedPosition{ positions->data() + data.startPositionIndex };
  const double* const queryPosition{ position.raw() };
  for (int i{ 0 }; i < totalPhaseSpaceDimension; ++i) {
    if (BumpDistribution::toDiscrete(storedPosition[i], inverseSegmentSpacing) != BumpDistribution::toDiscrete(queryPosition[i], inverseSegmentSpacing)) {
      return false;
    }
  }

  return true;
}

bool BumpDistribution::BumpDataEquality::operator()(const BumpData& first, const BumpData& second) const {
  const double* const firstPosition{ positions->data() + first.startPositionIndex };
  CoordinateSpan<const double> position(firstPosition, 1, totalPhaseSpaceDimension);
  return inSameSector(position, second);
}

size_t BumpDistribution::BumpDataHasher::getHashForPosition(CoordinateSpan<const double> position) const {
  if (position.totalSize() != static_cast<size_t>(totalPhaseSpaceDimension)) {
    return 0;
  }

  size_t hash{ 0 };
  const double* const positionData{ position.raw() };
  for (int i{ 0 }; i < totalPhaseSpaceDimension; ++i) {
    hash ^= std::hash<int>{}(BumpDistribution::toDiscrete(positionData[i], inverseSegmentSpacing)) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
  }
  return hash;
}

size_t BumpDistribution::BumpDataHasher::operator()(const BumpData& data) const {
  const double* const positionData{ positions->data() + data.startPositionIndex };
  CoordinateSpan<const double> position(positionData, 1, totalPhaseSpaceDimension);
  return getHashForPosition(position);
}
