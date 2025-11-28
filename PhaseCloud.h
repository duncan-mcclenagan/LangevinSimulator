#pragma once

#include <vector>

class PhaseCloud {

private:
  std::vector<double> data;
  int numberOfPoints;
  int numberOfParticles;
  int spatialDimension;
  int halfSize;

public:
  PhaseCloud(int numberOfPoints, int numberOfParticles, int spatialDimension);

  inline double& positionAt(size_t pointIndex, size_t particleIndex, size_t spatialCoordinateIndex) noexcept {
    return data[spatialCoordinateIndex * (numberOfParticles * numberOfPoints) +
      particleIndex * numberOfPoints +
      pointIndex];
  }

  inline double& velocityAt(size_t pointIndex, size_t particleIndex, size_t spatialCoordinateIndex) noexcept {
    return data[halfSize + (spatialCoordinateIndex * (numberOfParticles * numberOfPoints) +
      particleIndex * numberOfPoints +
      pointIndex)];
  }

  inline double* raw() noexcept {
    return data.data();
  }

};
