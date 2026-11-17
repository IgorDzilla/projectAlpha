#pragma once

#include "../models/EulerState.hpp"

// Два постоянных состояния по сторонам x0. Создаёт консервативные средние.
class RiemannProblem1D {
 public:
  RiemannProblem1D(EulerPrimitiveState left, EulerPrimitiveState right,
                   double discontinuity, double gamma = 1.4);
  State<1> initialState(const CartesianMesh<1>& mesh, double time = 0.0) const;
  State<1> initialState(const CartesianMesh<1>&&, double = 0.0) const = delete;

 private:
  EulerEquations gas_;
  EulerConservativeState left_, right_;
  double x0_;
};