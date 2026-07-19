#include <stdexcept>
#include "PhaseCloud.h"

PhaseCloud::PhaseCloud(int numberOfPoints, int spatialDimension, int numberOfParticles)
  : numberOfPoints(numberOfPoints), numberOfParticles(numberOfParticles), spatialDimension(spatialDimension), numCoordinates(numberOfPoints * spatialDimension * numberOfParticles), pointSize(spatialDimension * numberOfParticles)
{
  if (numberOfPoints <= 0) {
    throw std::invalid_argument("numberOfPoints must be positive");
  }
  if (spatialDimension <= 0) {
    throw std::invalid_argument("spatialDimension must be positive");
  }
  if (numberOfParticles <= 0) {
    throw std::invalid_argument("numberOfParticles must be positive");
  }

  data = std::vector<double>(numberOfPoints * spatialDimension * numberOfParticles * 2, 0);
}
