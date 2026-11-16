#include "numerics/methods/FiniteVolume.hpp"

#include <algorithm>
#include <sstream>

namespace {
std::string context(const char* what, std::size_t id, double time) {
  std::ostringstream out;
  out.precision(17);
  out << what << ' ' << id << " at t=" << time << ": ";
  return out.str();
}
}

FiniteVolume::FiniteVolume(const CartesianMesh<1>& mesh, double gamma,
                         double cfl, ExactRiemannOptions options)
    : mesh_(mesh), gas_(gamma), flux_(gamma, options), cfl_(cfl) {
  if (!std::isfinite(cfl_) || cfl_ <= 0.0 || cfl_ > 0.5)
    throw std::invalid_argument("CFL must be in (0, 0.5]");
}

void FiniteVolume::validate(const State<1>& state) const {
  validateEulerLayout(state);
  if (&state.mesh() != &mesh_)
    throw std::invalid_argument("FiniteVolume state belongs to another mesh");
  for (std::size_t i = 0; i < mesh_.cellCount(); ++i) {
    try { gas_.toPrimitive(eulerCell(state, i)); }
    catch (const std::exception& e) {
      throw std::runtime_error(context("cell", i, state.time()) + e.what());
    }
  }
}

void FiniteVolume::prepare(const State<1>& state) const {
  validate(state);
  const auto n = mesh_.cellCount();
  bool same = cachedCells_.size() == n && cachedTime_ == state.time();
  for (std::size_t i = 0; same && i < n; ++i) {
    const auto u = eulerCell(state, i);
    const auto& old = cachedCells_[i];
    same = u.density == old.density && u.momentum == old.momentum && u.energy == old.energy;
  }
  if (same) return;

  std::vector<EulerConservativeState> cells;
  std::vector<EulerPrimitiveState> primitive;
  cells.reserve(n); primitive.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    cells.push_back(eulerCell(state, i));
    primitive.push_back(gas_.toPrimitive(cells.back()));
  }
  std::vector<GodunovFaceResult> faces;
  faces.reserve(n + 1);
  double maxSpeed = 0.0;
  for (std::size_t face = 0; face <= n; ++face) {
    // На краях это пара одинаковых состояний — transmissive без хранимых ghost cells.
    const auto& left = primitive[face == 0 ? 0 : face - 1];
    const auto& right = primitive[face == n ? n - 1 : face];
    try { faces.push_back(flux_(left, right)); }
    catch (const std::exception& e) {
      throw std::runtime_error(context("face", face, state.time()) + e.what());
    }
    maxSpeed = std::max(maxSpeed, faces.back().signalSpeed);
  }
  if (!std::isfinite(maxSpeed) || maxSpeed <= 0.0)
    throw std::runtime_error("Invalid maximum wave speed");
  // Обновляем кэш только после успеха на всех гранях.
  cachedCells_ = std::move(cells);
  cachedFaces_ = std::move(faces);
  cachedTime_ = state.time();
  maxSpeed_ = maxSpeed;
}

double FiniteVolume::maxExplicitEulerStep(const State<1>& state) const {
  prepare(state);
  const double dt = cfl_ * (mesh_.step()[0] / maxSpeed_);
  if (!std::isfinite(dt) || dt <= 0.0)
    throw std::runtime_error("CFL time step cannot be represented");
  return dt;
}

void FiniteVolume::evaluate(const State<1>& state, State<1>::Fields& rates) const {
  prepare(state);
  if (rates.size() != 3)
    throw std::invalid_argument("Euler rates require three fields");
  for (const auto& rate : rates)
    if (&rate.mesh() != &mesh_ || rate.location() != FieldLocation::CellAverage ||
        rate.size() != mesh_.cellCount())
      throw std::invalid_argument("Euler rates have incompatible mesh or layout");
  for (std::size_t i = 0; i < mesh_.cellCount(); ++i) {
    const auto& a = cachedFaces_[i].flux;
    const auto& b = cachedFaces_[i+1].flux;
    // Одна грань используется обоими соседями с противоположными знаками.
    rates[0][i] = (a.mass - b.mass) / mesh_.step()[0];
    rates[1][i] = (a.momentum - b.momentum) / mesh_.step()[0];
    rates[2][i] = (a.energy - b.energy) / mesh_.step()[0];
  }
  state.validateFields(rates);
}
