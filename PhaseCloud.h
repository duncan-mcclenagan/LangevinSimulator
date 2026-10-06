#pragma once

#include <generator>
#include <vector>
#include "CoordinateSpan.h"

class PhaseCloud {

private:
  std::vector<double> data;

  inline int getPositionIndex(int pointIndex, int particleIndex, int spatialCoordinateIndex) const noexcept {
    return pointIndex * pointSize + particleIndex * spatialDimension + spatialCoordinateIndex;
  }

  inline int getVelocityIndex(int pointIndex, int particleIndex, int spatialCoordinateIndex) const noexcept {
    return numCoordinates + pointIndex * pointSize + particleIndex * spatialDimension + spatialCoordinateIndex;
  }

public:
  const int numCoordinates;
  const int numberOfPoints;
  const int numberOfParticles;
  const int spatialDimension;
  const int pointSize;

  PhaseCloud(int numberOfPoints, int numberOfParticles, int spatialDimension);

  struct PhasePoint {
    CoordinateSpan<double> position;
    CoordinateSpan<double> velocity;
  };

  std::generator<PhaseCloud::PhasePoint> enumeratePhasePoints() noexcept {
    double* position{ rawPosition() };
    double* velocity{ rawVelocity() };
    for (int i{ 0 }; i < numberOfPoints; ++i, position += pointSize, velocity += pointSize) {
      co_yield PhasePoint{
        CoordinateSpan<double>(position, numberOfParticles, spatialDimension),
        CoordinateSpan<double>(velocity, numberOfParticles, spatialDimension)
      };
    }
  }

  std::generator<CoordinateSpan<double>> enumeratePointPositions() noexcept {
    double* position{ rawPosition() };
    for (int i{ 0 }; i < numberOfPoints; ++i, position += pointSize) {
      co_yield CoordinateSpan<double>(position, numberOfParticles, spatialDimension);
    }
  }

  std::generator<CoordinateSpan<const double>> enumeratePointPositionsConst() const noexcept {
    const double* position{ rawPositionConst() };
    for (int i{ 0 }; i < numberOfPoints; ++i, position += pointSize) {
      co_yield CoordinateSpan<const double>(position, numberOfParticles, spatialDimension);
    }
  }

  std::generator<CoordinateSpan<double>> enumeratePointVelocities() noexcept {
    double* velocity{ rawVelocity() };
    for (int i{ 0 }; i < numberOfPoints; ++i, velocity += pointSize) {
      co_yield CoordinateSpan<double>(velocity, numberOfParticles, spatialDimension);
    }
  }

  std::generator<CoordinateSpan<const double>> enumeratePointVelocitiesConst() const noexcept {
    const double* velocity{ rawVelocityConst() };
    for (int i{ 0 }; i < numberOfPoints; ++i, velocity += pointSize) {
      co_yield CoordinateSpan<const double>(velocity, numberOfParticles, spatialDimension);
    }
  }

  inline CoordinateSpan<double> positionAt(int pointIndex) noexcept {
    return CoordinateSpan<double>(rawPosition(pointIndex), numberOfParticles, spatialDimension);
  }

  inline CoordinateSpan<const double> positionAtConst(int pointIndex) const noexcept {
    return CoordinateSpan<const double>(rawPositionConst(pointIndex), numberOfParticles, spatialDimension);
  }

  inline double& positionAt(int pointIndex, int particleIndex, int spatialCoordinateIndex) noexcept {
    return data[getPositionIndex(pointIndex, particleIndex, spatialCoordinateIndex)];
  }

  inline const double& positionAtConst(int pointIndex, int particleIndex, int spatialCoordinateIndex) const noexcept {
    return data[getPositionIndex(pointIndex, particleIndex, spatialCoordinateIndex)];
  }

  inline CoordinateSpan<double> velocityAt(int pointIndex) noexcept {
    return CoordinateSpan<double>(rawVelocity(pointIndex), numberOfParticles, spatialDimension);
  }

  inline CoordinateSpan<const double> velocityAtConst(int pointIndex) const noexcept {
    return CoordinateSpan<const double>(rawVelocityConst(pointIndex), numberOfParticles, spatialDimension);
  }

  inline double& velocityAt(int pointIndex, int particleIndex, int spatialCoordinateIndex) noexcept {
    return data[getVelocityIndex(pointIndex, particleIndex, spatialCoordinateIndex)];
  }

  inline const double& velocityAtConst(int pointIndex, int particleIndex, int spatialCoordinateIndex) const noexcept {
    return data[getVelocityIndex(pointIndex, particleIndex, spatialCoordinateIndex)];
  }

  inline double* rawPosition(int pointIndex = 0, int particleIndex = 0, int spatialCoordinateIndex = 0) noexcept {
    return data.data() + getPositionIndex(pointIndex, particleIndex, spatialCoordinateIndex);
  }

  inline const double* rawPositionConst(int pointIndex = 0, int particleIndex = 0, int spatialCoordinateIndex = 0) const noexcept {
    return data.data() + getPositionIndex(pointIndex, particleIndex, spatialCoordinateIndex);
  }

  inline double* rawVelocity(int pointIndex = 0, int particleIndex = 0, int spatialCoordinateIndex = 0) noexcept {
    return data.data() + getVelocityIndex(pointIndex, particleIndex, spatialCoordinateIndex);
  }

  inline const double* rawVelocityConst(int pointIndex = 0, int particleIndex = 0, int spatialCoordinateIndex = 0) const noexcept {
    return data.data() + getVelocityIndex(pointIndex, particleIndex, spatialCoordinateIndex);
  }

};
