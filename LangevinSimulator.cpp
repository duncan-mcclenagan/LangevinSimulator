#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
#include "BumpDistribution.h"
#include "EntropyCalculator.h"
#include "PhaseCloud.h"
#include "CoordinateSpan.h"

void runSimulation();
PhaseCloud createPhaseCloud(int numPoints, double velocitySpread, double positionSpread, int numParticles, int spatialDimension);
void computeForce(const BumpDistribution& bumps, CoordinateSpan<double> force, CoordinateSpan<double> position);

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

    int totalDimension{ static_cast<int>(position.totalSize()) };
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
  constexpr int numParticles{ 2 };
  
  constexpr double dampingCoefficient{ 0.03 };
  constexpr double beta{ 2E5 };

  constexpr double velocitySpread{ 0.5 };
  constexpr double positionSpread{ 1 };
  auto phaseCloud{ createPhaseCloud(numPoints, velocitySpread, positionSpread, numParticles, spatialDimension) };

  constexpr double segmentSpacing{ 0.1 };
  constexpr int numSegmentsInOneDimension{ 100 };
  constexpr double bumpMaxWidth{ 0.05 };
  constexpr double bumpFuzz{ 0.3 };
  constexpr double bumpHeight{ 1 };
  constexpr double bumpProbability{ 0.2 };
  auto bumps{ BumpDistribution::generate(numParticles,
    spatialDimension,
    segmentSpacing,
    numSegmentsInOneDimension,
    bumpProbability,
    BumpGenerator(segmentSpacing, bumpMaxWidth, bumpFuzz, bumpHeight)
  ) };

  EntropyCalculator entropyCalculator{};
  const double initialX1Entropy{ entropyCalculator.calculateEntropy(phaseCloud, bumps, 0, 0) };
  const double initialY2Entropy{ entropyCalculator.calculateEntropy(phaseCloud, bumps, 1, 1) };
  std::cout << "Initial particle 1 X entropy: " << initialX1Entropy << "\n";
  std::cout << "Initial particle 2 Y entropy: " << initialY2Entropy << "\n";

  std::random_device randomDevice;
  std::mt19937_64 generator(randomDevice());
  std::normal_distribution<double> normalDist(0, 1);
  std::vector<double> forceBuffer(numParticles * spatialDimension);
  CoordinateSpan<double> force(forceBuffer.data(), numParticles, spatialDimension);

  constexpr double halfTimeStep{ timeStep * 0.5 };
  const double firstVelocityDiffusionCoeff{ std::exp(-dampingCoefficient * halfTimeStep) };
  const double secondVelocityDiffusionCoeff{ std::sqrt((1 - firstVelocityDiffusionCoeff * firstVelocityDiffusionCoeff) / (beta)) };

  int totalPhaseSpaceDimension{ numParticles * spatialDimension };
  for (int i{ 0 }; i < numSteps; ++i) {
    // @todo Check this algorithm...
    for (auto [position, velocity] : phaseCloud.enumeratePhasePoints()) {
      // First reversible force update.
      computeForce(bumps, force, position);

      double* const rawPositionStart{ position.raw() };
      double* const rawVelocityStart{ velocity.raw() };
      double* const rawForceStart{ force.raw() };

      double* rawPosition{ rawPositionStart };
      double* rawVelocity{ rawVelocityStart };
      double* rawForce{ rawForceStart };
      for (int j{ 0 }; j < totalPhaseSpaceDimension; ++j, ++rawPosition, ++rawVelocity, ++rawForce) {
        // First diffusive velocity diffusion update.
        *rawVelocity = firstVelocityDiffusionCoeff * *rawVelocity + secondVelocityDiffusionCoeff * normalDist(generator);
        // First reversible position update.
        *rawPosition += halfTimeStep * *rawVelocity;
        // First reversible velocity update.
        *rawVelocity += halfTimeStep * *rawForce;
      }

      // Reversible force update. TODO: Should this be here?
      computeForce(bumps, force, position);

      rawPosition = rawPositionStart;
      rawVelocity = rawVelocityStart;
      rawForce = rawForceStart;
      for (int j{ 0 }; j < totalPhaseSpaceDimension; ++j, ++rawPosition, ++rawVelocity, ++rawForce) {
        // Second reversible position update.
        *rawPosition += halfTimeStep * *rawVelocity;
        // Second reversible velocity update.
        *rawVelocity += halfTimeStep * *rawForce;
        // Second diffusive velocity update.
        *rawVelocity = firstVelocityDiffusionCoeff * *rawVelocity + secondVelocityDiffusionCoeff * normalDist(generator);
      }
    }

    if (i % 100 == 0) {
      double currentXEntropy{ entropyCalculator.calculateEntropy(phaseCloud, bumps, 0, 0) };
      double currentYEntropy{ entropyCalculator.calculateEntropy(phaseCloud, bumps, 1, 1) };
      std::cout << "Time " << (i * timeStep) << ". Particle 1 X entropy: " << currentXEntropy << ". Particle 2 Y entropy: " << currentYEntropy << "\r";
    }
  }
}

void computeForce(const BumpDistribution& bumps, CoordinateSpan<double> force, CoordinateSpan<double> position) {
  bumps.computeForce(position.asConst(), force);
  // Add a constant "spring force" to the Y component of the second particle.
  force[1, 1] -= position[1, 1];
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
