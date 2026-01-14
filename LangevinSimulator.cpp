#include <cmath>
#include <iostream>
#include <random>
#include "PhaseCloud.h"

void runSimulation();

int main() {
  std::cout << "Welcome to the Langevin dynamics simulator!\n";
  runSimulation();
}

void runSimulation() {
  PhaseCloud phaseCloud(1, 1, 1);

  constexpr double timeStep{ 0.01 };
  constexpr int numSteps{ 10000 };
  
  constexpr double dampingCoefficient{ 0.01 };
  constexpr double beta{ 10 ^ 13 };

  constexpr double halfTimeStep{ timeStep * 0.5 };
  double firstVelocityDiffusionCoeff{ std::exp(-dampingCoefficient * halfTimeStep) };
  double secondVelocityDiffusionCoeff{ std::sqrt((1 - firstVelocityDiffusionCoeff * firstVelocityDiffusionCoeff) / (beta)) };

  std::random_device randomDevice;
  std::mt19937_64 generator(randomDevice());
  std::normal_distribution<double> normalDist(0, 1);

  auto& position{ phaseCloud.positionAt(0, 0, 0) };
  auto& velocity{ phaseCloud.velocityAt(0, 0, 0) };
  velocity = 1;
  position = 0;
  for (int i{ 0 }; i < numSteps; ++i) {
    const double force{ 0 };

    // First diffusive velocity diffusion update.
    velocity = firstVelocityDiffusionCoeff * velocity + secondVelocityDiffusionCoeff * normalDist(generator);
    // First reversible position update.
    position += halfTimeStep * velocity;
    // First reversible velocity update.
    velocity += halfTimeStep * force;
    // Reversible force update.
    // Second reversible position update.
    position += halfTimeStep * velocity;
    // Second reversible velocity update.
    velocity += halfTimeStep * force;
    // Second diffusive velocity update.
    velocity = firstVelocityDiffusionCoeff * velocity + secondVelocityDiffusionCoeff * normalDist(generator);
  }

  std::cout << "Final position: " << position << "\n";
  std::cout << "Final velocity: " << velocity << "\n";
}
