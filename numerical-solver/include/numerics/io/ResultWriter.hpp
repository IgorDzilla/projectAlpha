#pragma once

#include "../core/State.hpp"

// Формат вывода определит наследник. Solver подключается через Observer.
template <std::size_t Dim> class ResultWriter {
public:
  virtual ~ResultWriter() = default;
  virtual void write(const State<Dim> &state) = 0;
};
