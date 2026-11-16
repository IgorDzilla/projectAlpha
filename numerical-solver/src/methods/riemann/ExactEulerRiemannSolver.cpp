// Адаптировано из project Alpha/riemannSolver/src/exact-riemann-solver.cpp.
#include "numerics/methods/riemann/ExactEulerRiemannSolver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {

// Значение функции одной нелинейной волны f_K(p) и её производной.
struct WaveFunctionValue {
  double value;
  double derivative;
};

struct PressureGuess {
  double value;
  RiemannGuessType type;
};

double soundSpeed(const EulerPrimitiveState& state, double gamma) {
  return EulerEquations(gamma).soundSpeed(state);
}

// Функция f(p; W_K) имеет две ветви.
// p > p_K: ударная волна, используются соотношения Ранкина — Гюгонио.
// p <= p_K: волна разрежения, используются изэнтропические соотношения.
WaveFunctionValue evaluateWaveFunction(double pressure,
                                       const EulerPrimitiveState& state,
                                       double gamma) {
  if (!(pressure > 0.0) || !std::isfinite(pressure)) {
    throw std::domain_error("wave function requires positive pressure");
  }

  if (pressure > state.pressure) {
    const double A = 2.0 / ((gamma + 1.0) * state.density);
    const double B = (gamma - 1.0) / (gamma + 1.0) * state.pressure;
    const double root = std::sqrt(A / (pressure + B));

    return {
        (pressure - state.pressure) * root,
        root * (1.0 -
                (pressure - state.pressure) / (2.0 * (pressure + B))),
    };
  }

  const double a = soundSpeed(state, gamma);
  const double exponent = (gamma - 1.0) / (2.0 * gamma);
  const double ratio = pressure / state.pressure;

  return {
      2.0 * a / (gamma - 1.0) * (std::expm1(exponent * std::log(ratio))),
      std::pow(ratio, -(gamma + 1.0) / (2.0 * gamma)) /
          (state.density * a),
  };
}

// Правило GUESSP выбирает приближение, соответствующее ожидаемой конфигурации
// волн. Это уменьшает число итераций Ньютона на сильных разрывах.
PressureGuess guessPressure(const EulerPrimitiveState& left,
                            const EulerPrimitiveState& right,
                            double gamma) {
  const double aLeft = soundSpeed(left, gamma);
  const double aRight = soundSpeed(right, gamma);

  // Линеаризованная акустическая оценка (15.18).
  const double zLeft = left.density * aLeft;
  const double zRight = right.density * aRight;
  const double pLinear =
      (left.pressure * zRight + right.pressure * zLeft +
       (left.velocity - right.velocity) * zLeft * zRight) /
      (zLeft + zRight);

  const double pMin = std::min(left.pressure, right.pressure);
  const double pMax = std::max(left.pressure, right.pressure);

  if (pMax / pMin <= 2.0 && pLinear >= pMin && pLinear <= pMax) {
    return {pLinear, RiemannGuessType::Linearized};
  }

  if (pLinear < pMin) {
    // Предположение о двух волнах разрежения (15.19).
    const double z = (gamma - 1.0) / (2.0 * gamma);
    const double numerator =
        aLeft + aRight -
        0.5 * (gamma - 1.0) * (right.velocity - left.velocity);
    const double denominator =
        aLeft / std::pow(left.pressure, z) +
        aRight / std::pow(right.pressure, z);

    return {std::pow(numerator / denominator, 1.0 / z),
            RiemannGuessType::TwoRarefactions};
  }

  // Предположение о двух ударных волнах (15.20).
  // Давление pHat требуется только для приближённой оценки массовых потоков.
  const double pressureFloor =
      std::numeric_limits<double>::epsilon() * pMin;
  const double pHat = std::max(pLinear, pressureFloor);
  const double mLeft = std::sqrt(
      0.5 * left.density *
      ((gamma + 1.0) * pHat + (gamma - 1.0) * left.pressure));
  const double mRight = std::sqrt(
      0.5 * right.density *
      ((gamma + 1.0) * pHat + (gamma - 1.0) * right.pressure));

  const double estimate =
      (mLeft * mRight * (left.velocity - right.velocity) +
       mRight * left.pressure + mLeft * right.pressure) /
      (mLeft + mRight);

  return {std::max(estimate, pressureFloor), RiemannGuessType::TwoShocks};
}

double starDensity(const EulerPrimitiveState& state,
                   double pStar,
                   double gamma) {
  if (pStar > state.pressure) {
    // Плотность за ударной волной.
    const double ratio = pStar / state.pressure;
    const double g = (gamma - 1.0) / (gamma + 1.0);
    return state.density * (ratio + g) / (g * ratio + 1.0);
  }

  // Плотность за изэнтропической волной разрежения.
  return state.density *
         std::pow(pStar / state.pressure, 1.0 / gamma);
}

RiemannWaveInfo leftWaveInfo(const EulerPrimitiveState& left,
                      double pStar,
                      double uStar,
                      double gamma) {
  const double aLeft = soundSpeed(left, gamma);

  if (pStar > left.pressure) {
    const double speed =
        left.velocity -
        aLeft *
            std::sqrt((gamma + 1.0) / (2.0 * gamma) *
                          (pStar / left.pressure) +
                      (gamma - 1.0) / (2.0 * gamma));
    return {RiemannWaveType::Shock, speed, speed};
  }

  const double aStar =
      aLeft * std::pow(pStar / left.pressure,
                       (gamma - 1.0) / (2.0 * gamma));
  return {RiemannWaveType::Rarefaction,
          left.velocity - aLeft,  // внешняя граница веера
          uStar - aStar};         // граница веера у контактного разрыва
}

RiemannWaveInfo rightWaveInfo(const EulerPrimitiveState& right,
                       double pStar,
                       double uStar,
                       double gamma) {
  const double aRight = soundSpeed(right, gamma);

  if (pStar > right.pressure) {
    const double speed =
        right.velocity +
        aRight *
            std::sqrt((gamma + 1.0) / (2.0 * gamma) *
                          (pStar / right.pressure) +
                      (gamma - 1.0) / (2.0 * gamma));
    return {RiemannWaveType::Shock, speed, speed};
  }

  const double aStar =
      aRight * std::pow(pStar / right.pressure,
                        (gamma - 1.0) / (2.0 * gamma));
  return {RiemannWaveType::Rarefaction,
          right.velocity + aRight,  // внешняя граница веера
          uStar + aStar};           // граница веера у контактного разрыва
}

EulerPrimitiveState sampleLeftFan(const EulerPrimitiveState& left,
                             double S,
                             double gamma) {
  const double aLeft = soundSpeed(left, gamma);
  const double velocity =
      2.0 / (gamma + 1.0) *
      (aLeft + 0.5 * (gamma - 1.0) * left.velocity + S);
  const double a =
      2.0 / (gamma + 1.0) *
      (aLeft + 0.5 * (gamma - 1.0) * (left.velocity - S));
  const double soundRatio = a / aLeft;

  return {
      left.density * std::pow(soundRatio, 2.0 / (gamma - 1.0)),
      velocity,
      left.pressure *
          std::pow(soundRatio, 2.0 * gamma / (gamma - 1.0)),
  };
}

EulerPrimitiveState sampleRightFan(const EulerPrimitiveState& right,
                              double S,
                              double gamma) {
  const double aRight = soundSpeed(right, gamma);
  const double velocity =
      2.0 / (gamma + 1.0) *
      (-aRight + 0.5 * (gamma - 1.0) * right.velocity + S);
  const double a =
      2.0 / (gamma + 1.0) *
      (aRight - 0.5 * (gamma - 1.0) * (right.velocity - S));
  const double soundRatio = a / aRight;

  return {
      right.density * std::pow(soundRatio, 2.0 / (gamma - 1.0)),
      velocity,
      right.pressure *
          std::pow(soundRatio, 2.0 * gamma / (gamma - 1.0)),
  };
}

}  // namespace

ExactEulerRiemannSolver::ExactEulerRiemannSolver(double gamma,
                                               ExactRiemannOptions options)
    : gas_(gamma), options_(options) {
  if (!std::isfinite(options_.pressureRtol) || options_.pressureRtol <= 0.0 ||
      !std::isfinite(options_.residualRtol) || options_.residualRtol <= 0.0 ||
      options_.maxIterations == 0 || options_.maxBracketExpansions == 0) {
    throw std::invalid_argument("invalid exact Riemann options");
  }
}

RiemannSolution ExactEulerRiemannSolver::solve(
    const EulerPrimitiveState& left, const EulerPrimitiveState& right) const {
  gas_.validate(left);
  gas_.validate(right);
  const double gamma = gas_.gamma();
  const double aLeft = gas_.soundSpeed(left);
  const double aRight = gas_.soundSpeed(right);
  const double deltaVelocity = right.velocity - left.velocity;
  const double vacuumLimit = 2.0 * (aLeft + aRight) / (gamma - 1.0);
  if (!std::isfinite(deltaVelocity) || !std::isfinite(vacuumLimit)) {
    throw std::runtime_error("Riemann velocity scale overflow");
  }
  if (deltaVelocity >= vacuumLimit) {
    throw std::domain_error("vacuum is not supported by the exact Riemann solver");
  }
  const double velocityScale = std::max({aLeft, aRight, std::abs(deltaVelocity)});
  const auto residualAt = [&](double p) {
    const double residual = evaluateWaveFunction(p, left, gamma).value +
                            evaluateWaveFunction(p, right, gamma).value +
                            deltaVelocity;
    if (!std::isfinite(residual)) {
      throw std::runtime_error("non-finite Riemann pressure residual");
    }
    return residual;
  };

  RiemannSolution result;
  result.left = left;
  result.right = right;
  result.gamma = gamma;
  double pressure = left.pressure;
  // Равные состояния и чистый контакт решаются без итераций.
  if (left.pressure == right.pressure && left.velocity == right.velocity) {
    result.initialPressure = pressure;
    result.exactPressureRoot = true;
  } else {
    // g(0)<0 известно из проверки вакуума. Производную при p=0 не вычисляем.
    double low = 0.0;
    double high = std::max(left.pressure, right.pressure);
    double highResidual = residualAt(high);
    while (highResidual < 0.0) {
      if (result.bracketExpansions == options_.maxBracketExpansions) {
        throw std::runtime_error("Riemann pressure bracket expansion limit reached");
      }
      high *= 2.0;
      if (!std::isfinite(high)) {
        throw std::runtime_error("Riemann pressure bracket overflow");
      }
      ++result.bracketExpansions;
      highResidual = residualAt(high);
    }
    const PressureGuess guess = guessPressure(left, right, gamma);
    result.guessType = guess.type;
    pressure = (std::isfinite(guess.value) && guess.value > low && guess.value < high)
                   ? guess.value : low + 0.5 * (high - low);
    result.initialPressure = pressure;
    bool converged = false;
    if (highResidual == 0.0) {
      pressure = high;
      result.exactPressureRoot = true;
      converged = true;
    }
    for (std::size_t i = 0; !converged && i < options_.maxIterations; ++i) {
      const double residual = residualAt(pressure);
      result.iterations = i + 1;
      result.relativeResidual = std::abs(residual) / velocityScale;
      if (residual == 0.0) {
        result.relativeBracketWidth = 0.0;
        result.exactPressureRoot = true;
        converged = true;
        break;
      }
      if (residual < 0.0) low = pressure;
      else high = pressure;
      result.relativeBracketWidth = (high - low) / high;
      if (result.relativeBracketWidth <= options_.pressureRtol &&
          result.relativeResidual <= options_.residualRtol) {
        converged = true;
        break;
      }
      const double derivative = evaluateWaveFunction(pressure, left, gamma).derivative +
                                evaluateWaveFunction(pressure, right, gamma).derivative;
      double next = pressure - residual / derivative;
      // Ньютон принимается внутри центральных 80% интервала, иначе — середина.
      const double margin = 0.1 * (high - low);
      if (!std::isfinite(derivative) || derivative <= 0.0 || !std::isfinite(next) ||
          next <= low + margin || next >= high - margin) {
        next = low + 0.5 * (high - low);
      }
      if (!(next > low && next < high) || next == pressure) {
        throw std::runtime_error("Riemann pressure search stalled at floating-point precision");
      }
      pressure = next;
    }
    if (!converged) {
      throw std::runtime_error("Riemann pressure did not converge within maxIterations");
    }
  }

  const double fLeft = evaluateWaveFunction(pressure, left, gamma).value;
  const double fRight = evaluateWaveFunction(pressure, right, gamma).value;
  const double uStar = 0.5 * left.velocity + 0.5 * right.velocity +
                       0.5 * fRight - 0.5 * fLeft;
  result.pStar = pressure;
  result.uStar = uStar;
  result.starLeft = {starDensity(left, pressure, gamma), uStar, pressure};
  result.starRight = {starDensity(right, pressure, gamma), uStar, pressure};
  result.leftWave = leftWaveInfo(left, pressure, uStar, gamma);
  result.rightWave = rightWaveInfo(right, pressure, uStar, gamma);
  result.contactSpeed = uStar;
  gas_.validate(result.starLeft);
  gas_.validate(result.starRight);
  // Включаем и найденные фронты, и исходные акустические скорости.
  for (double speed : {result.leftWave.headSpeed, result.leftWave.tailSpeed,
                       result.rightWave.headSpeed, result.rightWave.tailSpeed,
                       uStar, std::abs(left.velocity) + aLeft,
                       std::abs(right.velocity) + aRight}) {
    if (!std::isfinite(speed)) {
      throw std::runtime_error("non-finite Riemann wave speed");
    }
    result.signalSpeed = std::max(result.signalSpeed, std::abs(speed));
  }
  return result;
}

EulerPrimitiveState ExactEulerRiemannSolver::sample(const RiemannSolution& solution,
                                          double S) const {
  if (!std::isfinite(S)) {
    throw std::invalid_argument("sample speed S must be finite");
  }
  if (!std::isfinite(solution.gamma) || std::abs(solution.gamma - gas_.gamma()) >
      16.0 * std::numeric_limits<double>::epsilon() * gas_.gamma()) {
    throw std::invalid_argument(
        "solution was produced with a different gamma");
  }

  gas_.validate(solution.left);
  gas_.validate(solution.right);
  gas_.validate(solution.starLeft);
  gas_.validate(solution.starRight);
  for (double speed : {solution.contactSpeed, solution.leftWave.headSpeed,
                       solution.leftWave.tailSpeed, solution.rightWave.headSpeed,
                       solution.rightWave.tailSpeed}) {
    if (!std::isfinite(speed)) {
      throw std::invalid_argument("invalid wave speed in Riemann solution");
    }
  }
  const auto checkedFan = [this](const EulerPrimitiveState& state) {
    gas_.validate(state);
    return state;
  };

  // Контакт движется со скоростью uStar. Сначала выбираем его сторону,
  // затем определяем положение луча относительно нелинейной волны.
  if (S < solution.contactSpeed) {
    if (solution.leftWave.type == RiemannWaveType::Shock) {
      return S < solution.leftWave.headSpeed ? solution.left
                                             : solution.starLeft;
    }

    if (S <= solution.leftWave.headSpeed) {
      return solution.left;
    }
    if (S >= solution.leftWave.tailSpeed) {
      return solution.starLeft;
    }
    return checkedFan(sampleLeftFan(solution.left, S, gas_.gamma()));
  }

  if (solution.rightWave.type == RiemannWaveType::Shock) {
    return S >= solution.rightWave.headSpeed ? solution.right
                                            : solution.starRight;
  }

  if (S >= solution.rightWave.headSpeed) {
    return solution.right;
  }
  if (S <= solution.rightWave.tailSpeed) {
    return solution.starRight;
  }
  return checkedFan(sampleRightFan(solution.right, S, gas_.gamma()));
}

std::string_view toString(RiemannWaveType type) {
  switch (type) {
    case RiemannWaveType::Shock:
      return "shock";
    case RiemannWaveType::Rarefaction:
      return "rarefaction";
  }
  return "unknown";
}

std::string_view toString(RiemannGuessType type) {
  switch (type) {
    case RiemannGuessType::Linearized:
      return "linearized";
    case RiemannGuessType::TwoRarefactions:
      return "two-rarefactions";
    case RiemannGuessType::TwoShocks:
      return "two-shocks";
  }
  return "unknown";
}

