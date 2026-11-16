#pragma once

#include "toml.hpp"
#include <cstdint>
#include <filesystem>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace tomlprocessor {

// Компилирует ограниченное математическое выражение один раз и затем
// вычисляет его для заданных координат. Копии разделяют скомпилированный код.
class MathExpression {
public:
  MathExpression(const std::string &formula,
                 const std::map<std::string, double> &parameters);

  // Ошибка синтаксиса выдаёт std::runtime_error при создании, а результат
  // вне области определения или бесконечность — при вычислении.
  double evaluate(double x, double t) const;

private:
  struct Impl;
  std::shared_ptr<Impl> impl_;
};

struct Problem {
  std::string type;
};

struct Parameters {
  std::map<std::string, double> values;
};

struct Domain {
  std::pair<double, double> x;
  std::pair<double, double> t;
};

struct Grid {
  double dx = 0.0;
  double dt = 0.0;
};

struct Solver {
  std::string method;
  bool check_stability = true;
  bool check_compatibility = true;
};

struct Segment {
  std::pair<double, double> interval;
  // Задаётся ровно одно из value и expression.
  std::optional<double> value;
  std::optional<std::string> expression;
  std::optional<MathExpression> compiled_expression;

  // Вычисляет значение сегмента; вызывающий выбирает сегмент по interval.
  double evaluate(double x, double t) const;
};

struct SourceFunction {
  std::string source; // "constant", "function", "table", "piecewise", "plugin"

  // Для constant.
  double value = 0.0;

  // Для function.
  std::string expression;
  std::optional<MathExpression> compiled_expression;

  // Для table. x_column — начальные условия, t_column — граничные.
  std::string path;
  std::string interpolation = "linear";
  std::string x_column;
  std::string t_column;
  std::string value_column;
  std::string delimiter = ",";
  std::string extrapolation = "error";

  // Для piecewise.
  std::vector<Segment> segments;

  // Для plugin.
  std::string module;
  std::string function;

  // Вычисляет только constant и function. Остальные источники требуют
  // отдельных загрузчиков и вызывают std::runtime_error.
  double evaluate(double x, double t) const;
};

struct Boundary : SourceFunction {
  std::string type; // "dirichlet", "neumann" или "newton"

  // Только для newton: alpha * u + beta * du/dn = g(t), beta != 0.
  // n — внешняя нормаль: du/dn = -u_x слева, +u_x справа.
  // g(t) задаётся полями SourceFunction. Оба коэффициента обязательны.
  std::optional<double> alpha;
  std::optional<double> beta;
};

using BoundariesLeft = Boundary;
using BoundariesRight = Boundary;
using InitialVelocity = SourceFunction;
using InitialPosition = SourceFunction;

struct Boundaries {
  BoundariesLeft left;
  BoundariesRight right;
};

struct Initials {
  InitialPosition position;
  InitialVelocity velocity;
};

struct Output {
  std::string directory = "results";
  std::string format = "csv";
  std::int64_t save_every = 1;
  std::vector<std::string> fields = {"position"};
};

struct Config {
  std::int64_t config_version = 1;
  Problem problem;
  Parameters parameters;
  Domain domain;
  Grid grid;
  Solver solver;
  Boundaries boundaries;
  Initials initials;
  std::optional<Output> output;
};

class ConfigProcessor {
public:
  // Чтение и структурная проверка TOML. Ошибки схемы: std::runtime_error;
  // ошибки синтаксиса TOML: toml::parse_error.
  // Выражения компилируются; CSV и плагины здесь не загружаются.
  static Config parse_file(const std::filesystem::path &path);
  static Config parse(const toml::table &table,
                      const std::filesystem::path &base_directory = ".");
};

} // namespace tomlprocessor
