#pragma once

// Одномерный идеальный газ. Эти типы общие для модели и локального решателя.
struct EulerPrimitiveState {
  double density = 0.0;  // rho — плотность.
  double velocity = 0.0; // u — скорость вдоль x.
  double pressure = 0.0; // p — давление.
};

struct EulerConservativeState {
  double density = 0.0;  // rho.
  double momentum = 0.0; // rho*u — импульс на единицу объёма.
  double energy = 0.0;   // E — полная энергия на единицу объёма, не на массу.
};

struct EulerFlux {
  double mass = 0.0;     // Поток массы: rho*u.
  double momentum = 0.0; // Поток импульса: rho*u*u + p.
  double energy = 0.0;   // Поток энергии: u*(E+p).
};

// Только физика одного состояния: без сетки, граничных условий и шага по времени.
class EulerEquations {
 public:
  explicit EulerEquations(double gamma = 1.4);
  double gamma() const noexcept { return gamma_; }

  void validate(const EulerPrimitiveState& state) const;
  EulerConservativeState toConservative(const EulerPrimitiveState& state) const;
  EulerPrimitiveState toPrimitive(const EulerConservativeState& state) const;
  double soundSpeed(const EulerPrimitiveState& state) const;
  EulerFlux physicalFlux(const EulerPrimitiveState& state) const;

 private:
  double gamma_;
};
