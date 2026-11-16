#include "numerics/methods/riemann/ExactEulerRiemannSolver.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {

struct ExpectedResult {
  double pStar;
  double uStar;
  double rhoStarLeft;
  double rhoStarRight;
  RiemannWaveType leftType;
  RiemannWaveType rightType;
};

struct ToroCase {
  const char* name;
  EulerPrimitiveState left;
  EulerPrimitiveState right;
  ExpectedResult expected;
};

// Табличные значения в задании округлены, поэтому сравниваем с небольшим
// абсолютным и относительным допуском, а не побитово.
bool nearlyEqual(double actual, double expected) {
  const double absoluteTolerance = 1.0e-5;
  const double relativeTolerance = 2.0e-5;
  return std::abs(actual - expected) <=
         absoluteTolerance +
             relativeTolerance * std::max(std::abs(actual), std::abs(expected));
}

bool checkValue(const std::string& label, double actual, double expected) {
  if (nearlyEqual(actual, expected)) {
    return true;
  }

  std::cerr << "    " << label << ": expected " << expected
            << ", got " << actual << '\n';
  return false;
}

bool runCase(const ToroCase& test, const ExactEulerRiemannSolver& solver) {
  const RiemannSolution result = solver.solve(test.left, test.right);

  bool passed = true;
  passed &= checkValue("p*", result.pStar, test.expected.pStar);
  passed &= checkValue("u*", result.uStar, test.expected.uStar);
  passed &= checkValue("rho*L", result.starLeft.density,
                       test.expected.rhoStarLeft);
  passed &= checkValue("rho*R", result.starRight.density,
                       test.expected.rhoStarRight);

  if (result.leftWave.type != test.expected.leftType) {
    std::cerr << "    left wave: expected "
              << toString(test.expected.leftType) << ", got "
              << toString(result.leftWave.type) << '\n';
    passed = false;
  }
  if (result.rightWave.type != test.expected.rightType) {
    std::cerr << "    right wave: expected "
              << toString(test.expected.rightType) << ", got "
              << toString(result.rightWave.type) << '\n';
    passed = false;
  }

  // Дополнительная проверка: отбор при S=0 и поток должны давать конечные
  // физически допустимые величины для каждого из пяти случаев.
  const EulerPrimitiveState face = solver.sample(result, 0.0);
  const EulerFlux flux = EulerEquations(solver.gamma()).physicalFlux(face);
  if (!(face.density > 0.0) || !(face.pressure > 0.0) ||
      !std::isfinite(flux.mass) || !std::isfinite(flux.momentum) ||
      !std::isfinite(flux.energy)) {
    std::cerr << "    sample(S=0) or flux is not physical\n";
    passed = false;
  }

  std::cout << (passed ? "[PASS] " : "[FAIL] ") << test.name
            << ", Pressure iterations: " << result.iterations
            << ", guess: " << toString(result.guessType) << '\n';
  return passed;
}

}  // namespace

int main() {
  const std::vector<ToroCase> cases{
      {
          "Toro 1 (Sod)",
          {1.0, 0.0, 1.0},
          {0.125, 0.0, 0.1},
          {0.30313, 0.92745, 0.42632, 0.26557,
           RiemannWaveType::Rarefaction, RiemannWaveType::Shock},
      },
      {
          "Toro 2 (two rarefactions)",
          {1.0, -2.0, 0.4},
          {1.0, 2.0, 0.4},
          {0.00189, 0.0, 0.02185, 0.02185,
           RiemannWaveType::Rarefaction,
           RiemannWaveType::Rarefaction},
      },
      {
          "Toro 3 (strong left pressure)",
          {1.0, 0.0, 1000.0},
          {1.0, 0.0, 0.01},
          {460.894, 19.5975, 0.57506, 5.99924,
           RiemannWaveType::Rarefaction, RiemannWaveType::Shock},
      },
      {
          "Toro 4 (strong right pressure)",
          {1.0, 0.0, 0.01},
          {1.0, 0.0, 100.0},
          {46.0950, -6.19633, 5.99242, 0.57511,
           RiemannWaveType::Shock, RiemannWaveType::Rarefaction},
      },
      {
          "Toro 5 (colliding shocks)",
          {5.99924, 19.5975, 460.894},
          {5.99242, -6.19633, 46.0950},
          {1691.64, 8.68975, 14.2823, 31.0426,
           RiemannWaveType::Shock, RiemannWaveType::Shock},
      },
  };

  try {
    const ExactEulerRiemannSolver solver(1.4);
    bool allPassed = true;
    for (const ToroCase& test : cases) {
      allPassed &= runCase(test, solver);
    }

    if (!allPassed) {
      std::cerr << "At least one Toro test failed.\n";
      return 1;
    }
    std::cout << "All five Toro tests passed.\n";
  } catch (const std::exception& error) {
    std::cerr << "Unexpected error: " << error.what() << '\n';
    return 1;
  }

  return 0;
}
