#pragma once

#include "../riemann/ExactEulerRiemannSolver.hpp"

struct GodunovFaceResult {
  EulerFlux flux;
  double signalSpeed; // Включает найденные фронты, а не только исходные |u|+a.
};

// Поток через одну неподвижную грань. Сетки и цикла по времени здесь нет.
class GodunovFlux {
 public:
  explicit GodunovFlux(double gamma = 1.4, ExactRiemannOptions options = {})
      : gas_(gamma), riemann_(gamma, options) {}

  GodunovFaceResult operator()(const EulerPrimitiveState& left,
                              const EulerPrimitiveState& right) const {
    const auto solution = riemann_.solve(left, right);
    return {gas_.physicalFlux(riemann_.sample(solution, 0.0)), solution.signalSpeed};
  }

 private:
  EulerEquations gas_;
  ExactEulerRiemannSolver riemann_;
};
