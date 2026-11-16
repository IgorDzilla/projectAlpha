#pragma once

#include "SpatialOperator.hpp"
#include "fluxes/GodunovFlux.hpp"
#include "../models/EulerState.hpp"

// Первый порядок для Эйлера в 1D: постоянное U в ячейке, Годунов на гранях.
// Обе границы transmissive: внешнее состояние равно ближайшему внутреннему.
// Сетка должна жить дольше оператора и не изменяться.
class FiniteVolume final : public SpatialOperator<1> {
 public:
  FiniteVolume(const CartesianMesh<1>& mesh, double gamma = 1.4,
               double cfl = 0.4, ExactRiemannOptions options = {});
  FiniteVolume(const CartesianMesh<1>&&, double = 1.4,
               double = 0.4, ExactRiemannOptions = {}) = delete;
  void validate(const State<1>& state) const override;
  double maxExplicitEulerStep(const State<1>& state) const override;
  void evaluate(const State<1>& state, State<1>::Fields& rates) const override;

 private:
  const CartesianMesh<1>& mesh_;
  EulerEquations gas_;
  GodunovFlux flux_;
  double cfl_;
  // Кэш сверяется со ВСЕМИ значениями U и временем, а не адресом State.
  mutable std::vector<EulerConservativeState> cachedCells_;
  mutable std::vector<GodunovFaceResult> cachedFaces_;
  mutable double cachedTime_ = 0.0;
  mutable double maxSpeed_ = 0.0;
  void prepare(const State<1>& state) const;
};
