#include "numerics/integrators/ExplicitEuler.hpp"
#include "numerics/io/ResultWriter.hpp"
#include "numerics/solvers/ExplicitSolver.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace {
void check(bool condition, const char *message) {
  if (!condition)
    throw std::runtime_error(message);
}

void close(double actual, double expected, const char *message) {
  check(std::abs(actual - expected) < 1e-12, message);
}

template <typename Error, typename Action>
void rejects(Action action, const char *message) {
  try {
    action();
  } catch (const Error &) {
    return;
  }
  throw std::runtime_error(message);
}

// Проверочная задача: dU_c/dt = (c+1)*(1+t) в каждой ячейке.
// Это проверка связки, а не реализация физической модели.
class Source final : public SpatialOperator<1> {
public:
  void evaluate(const State<1> &state, State<1>::Fields &rates) const override {
    for (std::size_t c = 0; c < rates.size(); ++c)
      rates[c].fill(static_cast<double>(c + 1) * (1.0 + state.time()));
  }
  double maxExplicitEulerStep(const State<1> &state) const override {
    return state.time() < 0.25 ? 0.25 : 0.125;
  }
};

class MemoryWriter final : public ResultWriter<1> {
public:
  std::vector<double> times;
  std::vector<double> values;
  void write(const State<1> &state) override {
    times.push_back(state.time());
    values.push_back(state.component(0)[0]);
  }
};

class FailingOperator final : public SpatialOperator<1> {
public:
  enum class Mode { Throw, Partial, WrongSize, Negative, NonFinite };
  explicit FailingOperator(Mode mode) : mode_(mode) {}

  void evaluate(const State<1> &, State<1>::Fields &rates) const override {
    if (mode_ == Mode::Throw)
      throw std::runtime_error("Operator failed");
    if (mode_ == Mode::WrongSize) {
      rates.clear();
      return;
    }
    if (mode_ == Mode::Partial) {
      rates[0][0] = 1.0;
      return;
    }
    for (auto &field : rates)
      field.fill(mode_ == Mode::Negative ? -100.0 :
                 std::numeric_limits<double>::infinity());
  }
  double maxExplicitEulerStep(const State<1> &) const override { return 1.0; }
  void validate(const State<1> &state) const override {
    if (state.component(0)[0] < 0.0)
      throw std::runtime_error("Physical state is invalid");
  }

private:
  Mode mode_;
};

class LimitOperator final : public SpatialOperator<1> {
  double limit_;

public:
  explicit LimitOperator(double limit) : limit_(limit) {}
  double maxExplicitEulerStep(const State<1> &) const override { return limit_; }
  void evaluate(const State<1> &, State<1>::Fields &rates) const override {
    for (auto &field : rates)
      field.fill(0.0);
  }
};

// Проверка, что Solver проверяет даже результат стороннего интегратора.
class BadIntegrator final : public TimeIntegrator<1> {
public:
  double maxTimeStep(const State<1> &, const SpatialOperator<1> &) const override {
    return 1.0;
  }
  State<1>::Fields advance(const State<1> &, const SpatialOperator<1> &,
                          double) const override {
    return {};
  }
};

void checkState() {
  CartesianMesh<1> mesh({0.0}, {0.5}, {2});
  State<1> state(mesh, {"density", "energy"}, 2.0);
  check(state.componentCount() == 2 && state.names()[1] == "energy", "State names");
  check(state.time() == 2.0 && state.component(0)[1] == 0.0, "State initialization");
  auto copy = state;
  copy.component(0)[0] = 3.0;
  check(state.component(0)[0] == 0.0, "Independent state copy");
  auto rates = state.makeFields(1.0);
  state.replaceFields(rates);
  check(state.component(1)[1] == 1.0 && state.time() == 2.0,
        "Replacing fields must preserve time");

  rejects<std::invalid_argument>([&] { State<1> bad(mesh, {}); }, "Empty components");
  rejects<std::invalid_argument>([&] { State<1> bad(mesh, {"u", "u"}); }, "Duplicate names");
  rejects<std::invalid_argument>([&] { State<1> bad(mesh, {""}); }, "Empty name");
  rejects<std::out_of_range>([&] { state.component(2); }, "Invalid component index");
  rejects<std::invalid_argument>([&] { state.replaceFields({}); }, "Invalid field count");
  check(state.component(1)[1] == 1.0, "Rejected fields must not change state");
  const double nan = std::numeric_limits<double>::quiet_NaN();
  rejects<std::invalid_argument>([&] { state.setTime(nan); }, "Non-finite time");
  check(state.time() == 2.0, "Rejected time must not change state");

  rates[0].fill(nan);
  rejects<std::runtime_error>([&] { state.replaceFields(rates); }, "Non-finite fields");
  CartesianMesh<1> otherMesh({0.0}, {0.5}, {2});
  copy.component(0) = Field<double, 1>(otherMesh);
  rejects<std::invalid_argument>([&] { copy.validate(); }, "Field rebound to wrong mesh");
  copy = state;
  copy.component(0) = Field<double, 1>(mesh, 0.0, FieldLocation::CellCenter);
  rejects<std::invalid_argument>([&] { copy.validate(); }, "Wrong field location");
}

void checkChain() {
  CartesianMesh<1> mesh({0.0}, {0.5}, {2});
  State<1> initial(mesh, {"a", "b"});
  ExplicitSolver<1> solver(initial, std::make_unique<Source>(),
                           std::make_unique<ExplicitEuler<1>>(), 0.3);
  Solver<1> &base = solver;
  MemoryWriter writer;
  base.run(0.6, [&](const State<1> &state) { writer.write(state); }, {3, 10});
  check(base.state().time() == 0.6, "Exact requested end time");
  close(base.state().component(0)[0], 0.728125, "Time-dependent source, component 0");
  close(base.state().component(1)[1], 1.45625, "Time-dependent source, component 1");
  check(initial.time() == 0.0 && initial.component(0)[0] == 0.0,
        "Solver must own independent initial values");
  check(writer.times.size() == 3 && writer.times[0] == 0.0 &&
        writer.times[1] == 0.5 && writer.times[2] == 0.6,
        "Initial, periodic and final output");
  close(writer.values.back(), base.state().component(0)[0], "Writer sees accepted state");

  base.run(0.725);
  close(base.state().component(0)[0], 0.928125, "Continuation from current state");
  rejects<std::invalid_argument>([&] { base.step(0.2); }, "Step above current limit");
  rejects<std::invalid_argument>([&] { base.run(0.0); }, "Backward run");
  for (double dt : {0.0, -0.1, std::numeric_limits<double>::infinity()})
    rejects<std::invalid_argument>([&] { base.step(dt); }, "Invalid step");

  MemoryWriter sameTime;
  base.run(base.state().time(), [&](const State<1> &s) { sameTime.write(s); });
  check(sameTime.times.size() == 1, "Zero-duration run outputs once");
}

void checkFailures() {
  CartesianMesh<1> mesh({0.0}, {0.5}, {2});
  State<1> initial(mesh, {"u"});
  initial.component(0).fill(1.0);
  for (auto mode : {FailingOperator::Mode::Throw, FailingOperator::Mode::Partial,
                    FailingOperator::Mode::WrongSize, FailingOperator::Mode::Negative,
                    FailingOperator::Mode::NonFinite}) {
    ExplicitSolver<1> solver(initial, std::make_unique<FailingOperator>(mode),
                             std::make_unique<ExplicitEuler<1>>(), 0.1);
    rejects<std::exception>([&] { solver.step(0.1); }, "Bad step must fail");
    check(solver.state().time() == 0.0 && solver.state().component(0)[0] == 1.0 &&
          solver.state().component(0)[1] == 1.0, "Rejected step changed accepted state");
  }
  ExplicitSolver<1> bad(initial, std::make_unique<Source>(),
                        std::make_unique<BadIntegrator>(), 0.1);
  rejects<std::invalid_argument>([&] { bad.step(0.1); }, "Invalid integrator output");
  check(bad.state().time() == 0.0, "Bad integrator changed time");

  for (double limit : {0.0, -1.0, std::numeric_limits<double>::quiet_NaN()}) {
    ExplicitSolver<1> solver(initial, std::make_unique<LimitOperator>(limit),
                             std::make_unique<ExplicitEuler<1>>(), 0.1);
    rejects<std::runtime_error>([&] { solver.maxTimeStep(); }, "Invalid stability limit");
  }
  ExplicitSolver<1> unrestricted(initial,
      std::make_unique<LimitOperator>(std::numeric_limits<double>::infinity()),
      std::make_unique<ExplicitEuler<1>>(), 0.1);
  check(unrestricted.maxTimeStep() == 0.1, "User limit must cap unrestricted operator");
  rejects<std::invalid_argument>([&] { unrestricted.run(1.0, {}, {0, 1}); },
                                 "Invalid output interval");
  rejects<std::invalid_argument>([&] { unrestricted.run(1.0, {}, {1, 0}); },
                                 "Invalid step budget");
  rejects<std::runtime_error>([&] { unrestricted.run(1.0, {}, {1, 1}); },
                              "Step budget must stop calculation");
  check(unrestricted.state().time() == 0.1, "Budget failure must preserve last accepted step");
  rejects<std::runtime_error>([&] {
    unrestricted.run(0.2, [](const State<1> &s) {
      if (s.time() == 0.2)
        throw std::runtime_error("Output failed");
    });
  }, "Output errors must propagate");
  check(unrestricted.state().time() == 0.2, "Output failure must not undo accepted step");

  State<1> late(mesh, {"u"}, 1e20);
  ExplicitSolver<1> stalled(late, std::make_unique<LimitOperator>(1.0),
                            std::make_unique<ExplicitEuler<1>>(), 1.0);
  rejects<std::runtime_error>([&] { stalled.step(0.1); }, "Unrepresentable time increment");
  check(stalled.state().time() == 1e20, "Stalled step changed time");
  rejects<std::invalid_argument>([&] {
    ExplicitSolver<1> missing(initial, nullptr, std::make_unique<ExplicitEuler<1>>(), 1.0);
  }, "Missing operator accepted");
  rejects<std::invalid_argument>([&] {
    ExplicitSolver<1> invalid(initial, std::make_unique<Source>(),
                              std::make_unique<ExplicitEuler<1>>(), 0.0);
  }, "Invalid user step limit accepted");
}

void checkNodalChain() {
  CartesianMesh<1> mesh({0.0}, {0.5}, {2});
  State<1> nodes(mesh, {"u"}, 0.0, FieldLocation::Node);
  check(nodes.component(0).size() == 3, "State must allocate all nodes");
  nodes.component(0)[2] = 4.0;
  auto rates = nodes.makeFields(2.0);
  check(rates[0].size() == 3 && rates[0].location() == FieldLocation::Node,
        "Nodal derivative buffers");
  nodes.validateFields(rates);
  auto wrong = rates;
  wrong[0] = Field<double, 1>(mesh, 0.0, FieldLocation::CellAverage);
  rejects<std::invalid_argument>([&] { nodes.validateFields(wrong); },
                                 "Cell fields accepted in nodal state");

  ExplicitSolver<1> solver(nodes, std::make_unique<Source>(),
                           std::make_unique<ExplicitEuler<1>>(), 0.25);
  solver.run(0.25);
  check(solver.state().location() == FieldLocation::Node &&
        solver.state().component(0).size() == 3, "Solver must preserve nodal layout");
  close(solver.state().component(0)[0], 0.25, "First node advanced");
  close(solver.state().component(0)[2], 4.25, "Last boundary node advanced");
  check(nodes.component(0)[2] == 4.0, "Original nodal state remains unchanged");

  CartesianMesh<2> plane({0.0, 0.0}, {1.0, 1.0}, {2, 3});
  State<2> nodes2D(plane, {"u", "v"}, 0.0, FieldLocation::Node);
  check(nodes2D.component(0).size() == 12 && nodes2D.component(1).size() == 12,
        "2D nodal state size");
  nodes2D.replaceFields(nodes2D.makeFields(1.0));
  nodes2D.validate();
}

static_assert(std::is_abstract_v<Solver<1>>, "Solver must remain abstract");
static_assert(std::has_virtual_destructor_v<Solver<1>>, "Solver needs virtual destruction");
static_assert(!std::is_constructible_v<State<1>, CartesianMesh<1> &&,
                                      std::vector<std::string>>, "Temporary mesh is unsafe");
} // namespace

int main() {
  try {
    checkState();
    checkChain();
    checkFailures();
    checkNodalChain();
    std::cout << "State and solver chain tests passed\n";
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
