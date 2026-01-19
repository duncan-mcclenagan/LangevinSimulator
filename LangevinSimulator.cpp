#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include "Bump.h"
#include "BumpDistribution.h"
#include "EntropyCalculator.h"
#include "PhaseCloud.h"

#include "Bumps.h"

void runSimulation();
double forceAt(double position);
PhaseCloud createPhaseCloud(int numPoints, double centerVelocity, double centerPosition, double velocitySpread, double positionSpread);

int main() {
  std::cout << "Welcome to the Langevin dynamics simulator!\n";
  runSimulation();
}

void runSimulation() {
  constexpr double timeStep{ 2E-4 };
  constexpr int numSteps{ 50000 };
  constexpr int numPoints{ 10000 };
  
  constexpr double dampingCoefficient{ 0.03 };
  constexpr double beta{ 2E5 };

  constexpr double halfTimeStep{ timeStep * 0.5 };
  const double firstVelocityDiffusionCoeff{ std::exp(-dampingCoefficient * halfTimeStep) };
  const double secondVelocityDiffusionCoeff{ std::sqrt((1 - firstVelocityDiffusionCoeff * firstVelocityDiffusionCoeff) / (beta)) };

  std::random_device randomDevice;
  std::mt19937_64 generator(randomDevice());
  std::normal_distribution<double> normalDist(0, 1);

  constexpr double centerVelocity{ 0.0 };
  constexpr double centerPosition{ 0.0 };
  constexpr double velocitySpread{ 0.5 };
  constexpr double positionSpread{ 1 };
  auto phaseCloud{ createPhaseCloud(numPoints, centerVelocity, centerPosition, velocitySpread, positionSpread) };
  EntropyCalculator entropyCalculator{};
  const double initialEntropy{ entropyCalculator.calculateEntropy(phaseCloud, bumps) };

  std::cout << "Initial entropy: " << initialEntropy << "\n";

  for (int i{ 0 }; i < numSteps; ++i) {
    double* position{ phaseCloud.rawPosition() };
    double* velocity{ phaseCloud.rawVelocity() };

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

  std::cout << "\nFinal entropy: " << entropyCalculator.calculateEntropy(phaseCloud, bumps) << "\n";
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

PhaseCloud createPhaseCloud(int numPoints, double centerVelocity, double centerPosition, double velocitySpread, double positionSpread) {
  std::random_device randomDevice;
  std::mt19937_64 generator(randomDevice());
  std::uniform_real_distribution<double> uniformDist(-1, 1);

  PhaseCloud phaseCloud(numPoints, 1, 1);
  double* position{ phaseCloud.rawPosition() };
  double* velocity{ phaseCloud.rawVelocity() };
  for (int i{ 0 }; i < numPoints; ++i, ++position, ++velocity) {
    *position = centerPosition + positionSpread * uniformDist(generator);
    *velocity = centerVelocity + velocitySpread * uniformDist(generator);
  }

  return phaseCloud;
}
