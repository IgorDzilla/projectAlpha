#include "toml-processor.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <initializer_list>
#include <set>
#include <stdexcept>
#include <string_view>

namespace tomlprocessor {
namespace {

using Keys = std::vector<std::string_view>;

[[noreturn]] void fail(const std::string &section, const std::string &message) {
  throw std::runtime_error("[" + section + "] " + message);
}

void check_keys(const toml::table &table, const std::string &section,
                const Keys &allowed) {
  for (const auto &[key, value] : table) {
    if (std::find(allowed.begin(), allowed.end(), key.str()) == allowed.end())
      fail(section, std::string(key.str()) + ": unknown or inapplicable field");
  }
}

const toml::node &required(const toml::table &table, std::string_view key,
                           const std::string &section) {
  const auto *node = table.get(key);
  if (!node)
    fail(section, std::string(key) + ": missing required field");
  return *node;
}

const toml::table &subtable(const toml::table &table, std::string_view key,
                            const std::string &section) {
  const auto *result = required(table, key, section).as_table();
  if (!result)
    fail(section, std::string(key) + ": expected table");
  return *result;
}

template <typename T>
T scalar(const toml::node &node, const std::string &section,
         std::string_view key) {
  auto value = node.value_exact<T>();
  if (!value)
    fail(section, std::string(key) + ": invalid value type");
  return *value;
}

template <typename T>
T get(const toml::table &table, std::string_view key,
      const std::string &section) {
  return scalar<T>(required(table, key, section), section, key);
}

template <typename T>
T get_or(const toml::table &table, std::string_view key,
         const std::string &section, T fallback) {
  const auto *node = table.get(key);
  return node ? scalar<T>(*node, section, key) : fallback;
}

std::string text(const toml::table &table, std::string_view key,
                 const std::string &section) {
  auto result = get<std::string>(table, key, section);
  if (result.find_first_not_of(" \t\r\n") == std::string::npos)
    fail(section, std::string(key) + ": must not be empty");
  return result;
}

double number(const toml::node &node, const std::string &section,
              std::string_view key) {
  if (!node.is_number())
    fail(section, std::string(key) + ": expected number");
  const double value = node.is_integer()
                           ? static_cast<double>(*node.value_exact<std::int64_t>())
                           : *node.value_exact<double>();
  if (!std::isfinite(value))
    fail(section, std::string(key) + ": must be finite");
  return value;
}

double number(const toml::table &table, std::string_view key,
              const std::string &section) {
  return number(required(table, key, section), section, key);
}

std::pair<double, double> interval(const toml::table &table,
                                   std::string_view key,
                                   const std::string &section) {
  const auto *array = required(table, key, section).as_array();
  if (!array || array->size() != 2)
    fail(section, std::string(key) + ": expected array of two numbers");
  const double a = number(*array->get(0), section, key);
  const double b = number(*array->get(1), section, key);
  if (b <= a)
    fail(section, std::string(key) + ": upper bound must exceed lower bound");
  return {a, b};
}

void choice(const std::string &value, std::initializer_list<std::string_view> choices,
            const std::string &section, std::string_view key) {
  if (std::find(choices.begin(), choices.end(), value) == choices.end())
    fail(section, std::string(key) + ": unsupported value \"" + value + "\"");
}

std::string resolve_path(const std::string &value,
                         const std::filesystem::path &base) {
  return std::filesystem::absolute(base / value).lexically_normal().string();
}

bool parameter_name(std::string_view name) {
  auto letter = [](char c) {
    return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || c == '_';
  };
  if (name.empty() || !letter(name.front()))
    return false;
  for (char c : name)
    if (!letter(c) && !(c >= '0' && c <= '9'))
      return false;
  // ExprTk сопоставляет имена без учёта регистра: X конфликтует с x.
  std::string folded(name);
  std::transform(folded.begin(), folded.end(), folded.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  const Keys reserved = {"x", "t", "pi", "e", "sin", "cos", "tan", "exp",
                         "log", "sqrt", "abs"};
  return std::find(reserved.begin(), reserved.end(), folded) == reserved.end();
}

// Разбирает один источник. Формулы компилируются здесь, поэтому неверное
// выражение получает диагностику с именем своей секции TOML.
SourceFunction read_source(const toml::table &table, const std::string &section,
                           bool boundary, const std::filesystem::path &base,
                           std::pair<double, double> range,
                           const std::map<std::string, double> &parameters) {
  SourceFunction result;
  result.source = text(table, "source", section);
  Keys keys = {"source"};
  if (boundary) {
    keys.push_back("type");
    if (get<std::string>(table, "type", section) == "newton") {
      keys.push_back("alpha");
      keys.push_back("beta");
    }
  }
  if (result.source == "constant") {
    keys.push_back("value");
    result.value = number(table, "value", section);
  } else if (result.source == "function") {
    keys.push_back("expression");
    result.expression = text(table, "expression", section);
    try {
      result.compiled_expression.emplace(result.expression, parameters);
    } catch (const std::runtime_error &error) {
      fail(section, "expression: " + std::string(error.what()));
    }
  } else if (result.source == "table") {
    keys.insert(keys.end(), {"path", "interpolation", "value_column", "delimiter",
                             "extrapolation", boundary ? "t_column" : "x_column"});
    result.path = resolve_path(text(table, "path", section), base);
    result.interpolation = get_or<std::string>(table, "interpolation", section, "linear");
    choice(result.interpolation, {"nearest", "linear", "cubic"}, section, "interpolation");
    result.value_column = text(table, "value_column", section);
    if (boundary)
      result.t_column = text(table, "t_column", section);
    else
      result.x_column = text(table, "x_column", section);
    result.delimiter = get_or<std::string>(table, "delimiter", section, ",");
    if (result.delimiter.size() != 1)
      fail(section, "delimiter: expected one character");
    result.extrapolation = get_or<std::string>(table, "extrapolation", section, "error");
    choice(result.extrapolation, {"error", "clamp"}, section, "extrapolation");
  } else if (result.source == "piecewise") {
    keys.push_back("segments");
    const auto *segments = required(table, "segments", section).as_array();
    if (!segments || segments->empty())
      fail(section, "segments: expected nonempty array of tables");
    for (std::size_t i = 0; i < segments->size(); ++i) {
      const std::string where = section + ".segments[" + std::to_string(i) + "]";
      const auto *segment = segments->get(i)->as_table();
      if (!segment)
        fail(where, "expected table");
      check_keys(*segment, where, {"interval", "value", "expression"});
      Segment part;
      part.interval = interval(*segment, "interval", where);
      if (segment->contains("value") == segment->contains("expression"))
        fail(where, "expected exactly one of value and expression");
      if (segment->contains("value"))
        part.value = number(*segment, "value", where);
      else {
        part.expression = text(*segment, "expression", where);
        try {
          part.compiled_expression.emplace(*part.expression, parameters);
        } catch (const std::runtime_error &error) {
          fail(where, "expression: " + std::string(error.what()));
        }
      }
      if (!result.segments.empty()) {
        const double previous_end = result.segments.back().interval.second;
        if (part.interval.first < previous_end)
          fail(where, "interval: segments overlap or are out of order");
        if (part.interval.first > previous_end && previous_end < range.second &&
            part.interval.first > range.first)
          fail(where, "interval: gap in the required range");
      }
      result.segments.push_back(std::move(part));
    }
    if (result.segments.front().interval.first > range.first ||
        result.segments.back().interval.second < range.second)
      fail(section, "segments: must cover the required range");
  } else if (result.source == "plugin") {
    keys.insert(keys.end(), {"module", "function"});
    result.module = text(table, "module", section);
    result.function = text(table, "function", section);
  } else {
    fail(section, "source: unsupported value \"" + result.source + "\"");
  }
  check_keys(table, section, keys);
  return result;
}

// Коэффициенты Ньютона проверяются отдельно от общего источника правой части.
Boundary read_boundary(const toml::table &table, const std::string &section,
                       const std::filesystem::path &base,
                       std::pair<double, double> range,
                       const std::map<std::string, double> &parameters) {
  Boundary result;
  result.type = text(table, "type", section);
  choice(result.type, {"dirichlet", "neumann", "newton"}, section, "type");
  if (result.type == "newton") {
    result.alpha = number(table, "alpha", section);
    result.beta = number(table, "beta", section);
    if (*result.beta == 0.0)
      fail(section, "beta: must be nonzero for newton");
  }
  static_cast<SourceFunction &>(result) = read_source(table, section, true, base, range, parameters);
  return result;
}

// Шаг должен делить длину области без изменения заданных пользователем границ.
void check_step(std::pair<double, double> range, double step, const char *key) {
  if (step <= 0.0)
    fail("grid", std::string(key) + ": must be > 0");
  const double count = (range.second - range.first) / step;
  // Допуск в числе интервалов, без изменения шага или границ области.
  if (!std::isfinite(count) || count < 1.0 ||
      std::abs(count - std::round(count)) > 1e-8)
    fail("grid", std::string(key) + ": domain length must be a multiple of step");
}

} // namespace

Config ConfigProcessor::parse_file(const std::filesystem::path &path) {
  // Относительные пути в конфигурации считаются от каталога самого файла.
  const auto absolute = std::filesystem::absolute(path);
  return parse(toml::parse_file(absolute.string()), absolute.parent_path());
}

Config ConfigProcessor::parse(const toml::table &table,
                              const std::filesystem::path &base_directory) {
  check_keys(table, "root", {"config_version", "problem", "parameters", "domain",
                             "grid", "solver", "boundaries", "initials", "output"});
  Config result;
  result.config_version = get<std::int64_t>(table, "config_version", "root");
  if (result.config_version != 1)
    fail("root", "config_version: only version 1 is supported");

  const auto &problem = subtable(table, "problem", "root");
  check_keys(problem, "problem", {"type"});
  result.problem.type = text(problem, "type", "problem");
  choice(result.problem.type, {"wave_1d"}, "problem", "type");

  const auto &parameters = subtable(table, "parameters", "root");
  std::set<std::string> folded_names;
  for (const auto &[key, node] : parameters) {
    if (!parameter_name(key.str()))
      fail("parameters", std::string(key.str()) + ": invalid or reserved name");
    std::string folded(key.str());
    std::transform(folded.begin(), folded.end(), folded.begin(), [](unsigned char c) {
      return static_cast<char>(std::tolower(c));
    });
    if (!folded_names.insert(folded).second)
      fail("parameters", std::string(key.str()) + ": name differs only by case");
    result.parameters.values.emplace(std::string(key.str()), number(node, "parameters", key.str()));
  }
  const double c = number(parameters, "c", "parameters");
  if (c <= 0.0)
    fail("parameters", "c: must be > 0");

  const auto &domain = subtable(table, "domain", "root");
  check_keys(domain, "domain", {"x", "t"});
  result.domain.x = interval(domain, "x", "domain");
  result.domain.t = interval(domain, "t", "domain");

  const auto &grid = subtable(table, "grid", "root");
  check_keys(grid, "grid", {"dx", "dt"});
  result.grid.dx = number(grid, "dx", "grid");
  result.grid.dt = number(grid, "dt", "grid");
  check_step(result.domain.x, result.grid.dx, "dx");
  check_step(result.domain.t, result.grid.dt, "dt");

  const auto &solver = subtable(table, "solver", "root");
  check_keys(solver, "solver", {"method", "check_stability", "check_compatibility"});
  result.solver.method = text(solver, "method", "solver");
  choice(result.solver.method, {"explicit_central"}, "solver", "method");
  result.solver.check_stability = get_or<bool>(solver, "check_stability", "solver", true);
  result.solver.check_compatibility = get_or<bool>(solver, "check_compatibility", "solver", true);
  if (result.solver.check_stability) {
    const long double cfl = static_cast<long double>(c) * result.grid.dt / result.grid.dx;
    if (!std::isfinite(cfl) || cfl > 1.0L)
      fail("solver", "explicit_central: CFL exceeds 1.0");
  }

  const auto &boundaries = subtable(table, "boundaries", "root");
  check_keys(boundaries, "boundaries", {"left", "right"});
  result.boundaries.left = read_boundary(subtable(boundaries, "left", "boundaries"),
      "boundaries.left", base_directory, result.domain.t, result.parameters.values);
  result.boundaries.right = read_boundary(subtable(boundaries, "right", "boundaries"),
      "boundaries.right", base_directory, result.domain.t, result.parameters.values);

  const auto &initials = subtable(table, "initials", "root");
  check_keys(initials, "initials", {"position", "velocity"});
  result.initials.position = read_source(subtable(initials, "position", "initials"),
      "initials.position", false, base_directory, result.domain.x, result.parameters.values);
  result.initials.velocity = read_source(subtable(initials, "velocity", "initials"),
      "initials.velocity", false, base_directory, result.domain.x, result.parameters.values);

  if (table.contains("output")) {
    const auto &output = subtable(table, "output", "root");
    check_keys(output, "output", {"directory", "format", "save_every", "fields"});
    Output settings;
    settings.directory = resolve_path(output.contains("directory")
        ? text(output, "directory", "output") : "results", base_directory);
    settings.format = get_or<std::string>(output, "format", "output", "csv");
    choice(settings.format, {"csv"}, "output", "format");
    settings.save_every = get_or<std::int64_t>(output, "save_every", "output", 1);
    if (settings.save_every < 1)
      fail("output", "save_every: must be >= 1");
    if (output.contains("fields")) {
      const auto *fields = output.get("fields")->as_array();
      if (!fields || fields->empty())
        fail("output", "fields: expected nonempty array of strings");
      settings.fields.clear();
      for (const auto &node : *fields) {
        const auto field = scalar<std::string>(node, "output", "fields");
        choice(field, {"position", "velocity"}, "output", "fields");
        if (std::find(settings.fields.begin(), settings.fields.end(), field) != settings.fields.end())
          fail("output", "fields: duplicate field \"" + field + "\"");
        settings.fields.push_back(field);
      }
    }
    result.output = std::move(settings);
  }
  return result;
}

} // namespace tomlprocessor
