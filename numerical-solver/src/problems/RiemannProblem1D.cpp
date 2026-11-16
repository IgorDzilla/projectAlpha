#include "numerics/problems/RiemannProblem1D.hpp"

#include <algorithm>

RiemannProblem1D::RiemannProblem1D(EulerPrimitiveState left,
    EulerPrimitiveState right, double discontinuity, double gamma)
    : gas_(gamma), left_(gas_.toConservative(left)),
      right_(gas_.toConservative(right)), x0_(discontinuity) {
  if (!std::isfinite(x0_))
    throw std::invalid_argument("Riemann discontinuity must be finite");
}

State<1> RiemannProblem1D::initialState(const CartesianMesh<1>& mesh, double time) const {
  if (!(x0_ > mesh.origin()[0] && x0_ < mesh.upperBound()[0]))
    throw std::invalid_argument("Riemann discontinuity must lie inside the domain");
  State<1> result(mesh, {"density", "momentum", "energy"}, time);
  for (std::size_t i = 0; i < mesh.cellCount(); ++i) {
    const double a = mesh.faceCoordinate(0, i);
    const double b = mesh.faceCoordinate(0, i + 1);
    if (!(b > a))
      throw std::invalid_argument("Mesh faces are indistinguishable at this precision");
    // При пересечении ячейки разрывом смешиваем U, а не rho/u/p по отдельности.
    const double fraction = std::clamp((x0_ - a) / (b - a), 0.0, 1.0);
    result.component(0)[i] = fraction * left_.density + (1-fraction) * right_.density;
    result.component(1)[i] = fraction * left_.momentum + (1-fraction) * right_.momentum;
    result.component(2)[i] = fraction * left_.energy + (1-fraction) * right_.energy;
    gas_.toPrimitive(eulerCell(result, i));
  }
  result.validate();
  return result;
}
