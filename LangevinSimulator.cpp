#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
#include "Bump.h"
#include "BumpDistribution.h"
#include "EntropyCalculator.h"
#include "PhaseCloud.h"

void runSimulation();
double forceAt(double position);
PhaseCloud createPhaseCloud(int numPoints, double velocitySpread, double positionSpread);

class BumpGenerator {

private:
  const double segmentHalfSpacing;
  const double bumpMaxWidth;
  double bumpPositionVariation;
  double bumpWidthVariationMultiplier;
  const double bumpHeight;

  std::mt19937_64 generator;
  std::uniform_real_distribution<double> uniformDist;

public:
  BumpGenerator(double segmentSpacing, double bumpMaxWidth, double bumpFuzz, double bumpHeight)
    : segmentHalfSpacing{ std::max(0.0, segmentSpacing * 0.5) },
    bumpMaxWidth{ std::min(bumpMaxWidth, segmentSpacing) },
    bumpHeight{ std::max(0.0, bumpHeight) },
    generator((std::random_device())()),
    uniformDist(-1, 1)
  {
    bumpFuzz = std::max(0.0, std::min(bumpFuzz, 1.0));
    bumpPositionVariation = segmentHalfSpacing * bumpFuzz;
    bumpWidthVariationMultiplier = bumpFuzz;
  }

  void operator()(double& width, double& height, CoordinateSpan<double> position) {
    width = bumpMaxWidth * bumpWidthVariationMultiplier * uniformDist(generator);
    height = bumpHeight;

    int totalDimension{ position.totalSize() };
    double* positionData{ position.raw() };
    for (int i =0; i < totalDimension; ++i) {
      positionData[i] = bumpPositionVariation * uniformDist(generator);
    }
  }

};

int main() {
  std::cout << "Welcome to the Langevin dynamics simulator!\n";
  runSimulation();

  return 0;
}

void runSimulation() {
  constexpr double timeStep{ 2E-4 };
  constexpr int numSteps{ 50000 };
  constexpr int numPoints{ 10000 };
  constexpr int spatialDimension{ 2 };
  constexpr int numParticles{ 1 };
  
  constexpr double dampingCoefficient{ 0.03 };
  constexpr double beta{ 2E5 };

  constexpr double halfTimeStep{ timeStep * 0.5 };
  const double firstVelocityDiffusionCoeff{ std::exp(-dampingCoefficient * halfTimeStep) };
  const double secondVelocityDiffusionCoeff{ std::sqrt((1 - firstVelocityDiffusionCoeff * firstVelocityDiffusionCoeff) / (beta)) };

  constexpr double velocitySpread{ 0.5 };
  constexpr double positionSpread{ 1 };
  auto phaseCloud{ createPhaseCloud(numPoints, velocitySpread, positionSpread) };

  constexpr double segmentSpacing{ 0.1 };
  constexpr int numSegmentsInOneDimension{ 100 };
  constexpr double bumpMaxWidth{ 0.05 };
  constexpr double bumpFuzz{ 0.3 };
  constexpr double bumpHeight{ 1 };
  auto bumps{ BumpDistribution::generate(numParticles, spatialDimension, segmentSpacing, numSegmentsInOneDimension, BumpGenerator(segmentSpacing, bumpMaxWidth, bumpFuzz, bumpHeight)) };

  EntropyCalculator entropyCalculator{};
  const double initialEntropy{ entropyCalculator.calculateEntropy(phaseCloud, bumps) };
  std::cout << "Initial entropy: " << initialEntropy << "\n";

  std::random_device randomDevice;
  std::mt19937_64 generator(randomDevice());
  std::normal_distribution<double> normalDist(0, 1);
  for (int i{ 0 }; i < numSteps; ++i) {
    // @todo Check this algorithm...
    for (auto [position, velocity] : phaseCloud.enumeratePhasePoints()) {

    }
    for (int j{ 0 }; j < phaseCloud.numberOfPoints; ++j, ++position, ++velocity) {
      double force{ forceAt(*position) };

      // First diffusive velocity diffusion update.
      *velocity = firstVelocityDiffusionCoeff * *velocity + secondVelocityDiffusionCoeff * normalDist(generator);
      // First reversible position update.
      *position += halfTimeStep * *velocity;
      // First reversible velocity update.
      *velocity += halfTimeStep * force;
      // Reversible force update. TODO: Should this be here?
      force = forceAt(*position);
      // Second reversible position update.
      *position += halfTimeStep * *velocity;
      // Second reversible velocity update.
      *velocity += halfTimeStep * force;
      // Second diffusive velocity update.
      *velocity = firstVelocityDiffusionCoeff * *velocity + secondVelocityDiffusionCoeff * normalDist(generator);
    }

    if (i % 100 == 0) {
      double currentEntropy{ entropyCalculator.calculateEntropy(phaseCloud, bumps) };
      std::cout << "Time " << (i * timeStep) << ". Entropy: " << currentEntropy << "\r";
    }
  }
}

double forceAt(double position) {
  // Use the bump function V(x) = -H*e^(1/(((x - a)/w)^2 - 1)) to produce our minimum.
  // H is the height of the bump, w is the width, and a is the center.
  // The corresponding force is -H*e^(1/(((x - a)/w)^2 - 1))*(2*(x - a)/w^2)/(((x - a)/w)^2 - 1)^2.
  // We divide a region up into segments and place a bump in each segment with some (randomzied but preset)
  // center offset and width (we keep the height constant for now).

  constexpr double height{ 0.2 };
  constexpr double segmentSpacing{ 0.1 };
  constexpr double inverseSegmentSpacing{ 1 / segmentSpacing };
  constexpr int numSegments{ 201 };
  constexpr int centerSegmentIndex{ numSegments / 2 };

  const auto bump{ bumps.getNearestBump(position) };
  const double relativePosition{ position - bump.center };

  // See if we are within the bump.
  if (std::abs(relativePosition) > bump.width) {
    return 0;
  }

  // Calculate the force using the formula above.
  // -H*e^(1/(((x - a)/w)^2 - 1))*(2*(x - a)/w^2)/(((x - a)/w)^2 - 1)^2
  const double inverseWidth{ 1 / bump.width };
  const double normalizedPosition{ relativePosition * inverseWidth };
  const double expArg{ 1 / (normalizedPosition * normalizedPosition - 1) };

  return -height * std::exp(expArg) * (2 * relativePosition * inverseWidth * inverseWidth) * expArg * expArg;
}

PhaseCloud createPhaseCloud(int numPoints, double velocitySpread, double positionSpread, int numParticles, int spatialDimension) {
  std::random_device randomDevice;
  std::mt19937_64 generator(randomDevice());
  std::uniform_real_distribution<double> uniformDist(-1, 1);

  PhaseCloud phaseCloud(numPoints, numParticles, spatialDimension);
  double* position{ phaseCloud.rawPosition() };
  double* velocity{ phaseCloud.rawVelocity() };
  for (int i{ 0 }; i < phaseCloud.numCoordinates; ++i, ++position, ++velocity) {
    *position = positionSpread * uniformDist(generator);
    *velocity = velocitySpread * uniformDist(generator);
  }

  return phaseCloud;
}
