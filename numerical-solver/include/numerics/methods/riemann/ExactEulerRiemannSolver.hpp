#pragma once

#include "numerics/models/EulerEquations.hpp"

#include <cstddef>
#include <string_view>

enum class RiemannWaveType {
  Shock,       // Ударная волна: скачок состояния.
  Rarefaction, // Разрежение: непрерывный веер; также обозначает волну нулевой силы.
};

// Как выбрано первое приближение давления; это диагностика, не тип решателя.
enum class RiemannGuessType {
  Linearized,
  TwoRarefactions,
  TwoShocks,
};

// Для волны разрежения headSpeed и tailSpeed различны:
// headSpeed — внешняя граница веера, tailSpeed — граница у звёздной области.
// Для ударной волны обе скорости равны скорости разрыва.
struct RiemannWaveInfo {
  RiemannWaveType type = RiemannWaveType::Rarefaction;
  double headSpeed = 0.0;
  double tailSpeed = 0.0;
};

// Полный результат решения одной задачи Римана.
// В него включены исходные состояния, чтобы sample() мог затем отбирать
// решение на любом луче S = x/t без повторного решения уравнения для pStar.
struct RiemannSolution {
  EulerPrimitiveState left;
  EulerPrimitiveState right;
  EulerPrimitiveState starLeft;  // Промежуточное состояние слева от контакта.
  EulerPrimitiveState starRight; // Промежуточное состояние справа от контакта.

  double gamma = 1.4;
  double pStar = 0.0; // Общее давление между внешними волнами.
  double uStar = 0.0; // Общая скорость газа по обе стороны контакта.

  RiemannWaveInfo leftWave;
  RiemannWaveInfo rightWave;
  double contactSpeed = 0.0;

  RiemannGuessType guessType = RiemannGuessType::Linearized;
  double initialPressure = 0.0;
  std::size_t iterations = 0; // Уточнения давления: Ньютон и деление интервала.
  std::size_t bracketExpansions = 0;
  double relativeResidual = 0.0; // Невязка уравнения давления / масштаб скорости.
  double relativeBracketWidth = 0.0;
  bool exactPressureRoot = false; // Невязка точно нулевая в арифметике double.
  double signalSpeed = 0.0; // Максимальный модуль скорости для будущего CFL.
};

struct ExactRiemannOptions {
  double pressureRtol = 1.0e-10; // Относительная ширина интервала давления.
  double residualRtol = 1.0e-10; // Допустимая нормированная невязка.
  std::size_t maxIterations = 256;
  std::size_t maxBracketExpansions = 128;
};

// Локальная задача Эйлера в 1D. Не является Solver<1>: не хранит поле и время.
// Основа перенесена из пользовательского riemannSolver; алгоритм уточнения усилен.
class ExactEulerRiemannSolver {
 public:
  explicit ExactEulerRiemannSolver(double gamma = 1.4,
                                 ExactRiemannOptions options = {});

  // W слева/справа от одной грани. Результат владеет копиями входных данных.
  // Вакуум и несходимость дают исключение, подмены приближённым решателем нет.
  RiemannSolution solve(const EulerPrimitiveState& left,
                        const EulerPrimitiveState& right) const;

  // xi = (x-x_face)/(t-t_n); для неподвижной грани xi=0.
  // Принимать результат solve() с тем же gamma; поля результата не менять.
  // На скачке возвращается правый пространственный след.
  EulerPrimitiveState sample(const RiemannSolution& solution, double xi) const;

  double gamma() const noexcept { return gas_.gamma(); }

 private:
  EulerEquations gas_;
  ExactRiemannOptions options_;
};

std::string_view toString(RiemannWaveType type);
std::string_view toString(RiemannGuessType type);
