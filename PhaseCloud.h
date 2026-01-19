#pragma once

#include <vector>

class PhaseCloud {

private:
  std::vector<double> data;
  int halfSize;

  inline int getPositionIndex(int pointIndex, int particleIndex, int spatialCoordinateIndex) const noexcept {
    return spatialCoordinateIndex * (numberOfParticles * numberOfPoints) + particleIndex * numberOfPoints + pointIndex;
  }

  inline int getVelocityIndex(int pointIndex, int particleIndex, int spatialCoordinateIndex) const noexcept {
    return halfSize + (spatialCoordinateIndex * (numberOfParticles * numberOfPoints) + particleIndex * numberOfPoints + pointIndex);
  }

public:
  const int numberOfPoints;
  const int numberOfParticles;
  const int spatialDimension;

  PhaseCloud(int numberOfPoints, int numberOfParticles, int spatialDimension);

  // Rely on automatic move constructor.

  inline double& positionAt(int pointIndex, int particleIndex, int spatialCoordinateIndex) noexcept {
    return data[getPositionIndex(pointIndex, particleIndex, spatialCoordinateIndex)];
  }

  inline double& velocityAt(int pointIndex, int particleIndex, int spatialCoordinateIndex) noexcept {
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
