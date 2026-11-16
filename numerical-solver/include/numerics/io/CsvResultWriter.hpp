#pragma once

#include "ResultWriter.hpp"
#include "../models/EulerState.hpp"
#include <filesystem>
#include <fstream>

// Один solution.csv: строка на ячейку и сохранённый момент времени.
class CsvResultWriter final : public ResultWriter<1> {
 public:
  CsvResultWriter(const std::filesystem::path& directory, double gamma,
                  std::vector<std::string> fields);
  void write(const State<1>& state) override;
  std::size_t snapshots() const noexcept { return snapshots_; }
 private:
  EulerEquations gas_;
  std::vector<std::string> fields_;
  std::ofstream stream_;
  std::size_t snapshots_ = 0;
};
