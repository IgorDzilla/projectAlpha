#include "toml-processor.hpp"

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>

namespace {
const std::string input = R"toml(
config_version = 1
[problem]
type = "wave_1d"
[parameters]
c = 1
L = 10.0
[domain]
x = [0, 10]
t = [0, 20]
[grid]
dx = 0.1
dt = 0.05
[solver]
method = "explicit_central"
[boundaries.left]
type = "dirichlet"
source = "constant"
value = 0
[boundaries.right]
type = "newton"
alpha = 2
beta = 0.5
source = "function"
expression = "sin(t)"
[initials.position]
source = "piecewise"
[[initials.position.segments]]
interval = [0, 5]
value = 0
[[initials.position.segments]]
interval = [5, 10]
expression = "x - 5"
[initials.velocity]
source = "table"
path = "data/velocity.csv"
x_column = "x"
value_column = "velocity"
[output]
)toml";

void check(bool condition, const std::string &message) {
  if (!condition)
    throw std::runtime_error(message);
}

std::string replace(std::string source, const std::string &from, const std::string &to) {
  const auto pos = source.find(from);
  check(pos != std::string::npos, "test replacement not found: " + from);
  source.replace(pos, from.size(), to);
  return source;
}

void rejects(const std::string &source, const std::string &diagnostic) {
  try {
    tomlprocessor::ConfigProcessor::parse(toml::parse(source));
  } catch (const std::runtime_error &error) {
    check(std::string(error.what()).find(diagnostic) != std::string::npos,
          "unexpected error: " + std::string(error.what()));
    return;
  }
  throw std::runtime_error("expected rejection: " + diagnostic);
}
} // namespace

int main() {
  using tomlprocessor::ConfigProcessor;
  const auto base = std::filesystem::temp_directory_path() / "toml-parser-base";
  const auto config = ConfigProcessor::parse(toml::parse(input), base);
  check(config.config_version == 1 && config.parameters.values.at("c") == 1,
        "version and numeric parameters");
  check(config.domain.x.second == 10 && config.grid.dt == 0.05, "domain and grid");
  check(config.solver.check_stability && config.solver.check_compatibility, "solver defaults");
  check(config.boundaries.right.alpha == 2 && config.boundaries.right.beta == 0.5,
        "newton coefficients");
  check(config.boundaries.right.expression == "sin(t)", "expression preserved");
  check(std::abs(config.boundaries.right.evaluate(0, std::acos(-1.0) / 2) - 1.0) < 1e-12,
        "boundary expression evaluates with t");
  const auto &segments = config.initials.position.segments;
  check(segments.size() == 2 && segments[0].value == 0 && !segments[0].expression &&
        !segments[1].value && segments[1].expression == "x - 5", "piecewise variants");
  check(segments[0].evaluate(2, 0) == 0 && segments[1].evaluate(7, 0) == 2,
        "piecewise expressions evaluate with x");
  const auto powered = replace(input, "expression = \"x - 5\"",
                               "expression = \"2**3 + L - x\"");
  check(ConfigProcessor::parse(toml::parse(powered)).initials.position.segments[1]
            .evaluate(7, 0) == 11, "power and parameter expression");
  const auto constants = replace(input, "expression = \"sin(t)\"",
                                 "expression = \"cos(pi) + exp(0) + log(e)\"");
  check(std::abs(ConfigProcessor::parse(toml::parse(constants))
                     .boundaries.right.evaluate(0, 0) - 1.0) < 1e-12,
        "allowed functions and constants");
  check(config.initials.velocity.path == (base / "data/velocity.csv").string(), "CSV path");
  check(config.initials.velocity.interpolation == "linear" &&
        config.initials.velocity.delimiter == "," &&
        config.initials.velocity.extrapolation == "error", "table defaults");
  check(config.output && config.output->directory == (base / "results").string() &&
        config.output->format == "csv" && config.output->save_every == 1 &&
        config.output->fields == std::vector<std::string>{"position"}, "output defaults");
  check(!ConfigProcessor::parse(toml::parse(replace(input, "[output]", ""))).output,
        "optional output");

  const auto plugin = replace(input, "source = \"function\"\nexpression = \"sin(t)\"",
      "source = \"plugin\"\nmodule = \"custom.data\"\nfunction = \"boundary\"");
  check(ConfigProcessor::parse(toml::parse(plugin)).boundaries.right.module == "custom.data",
        "plugin metadata");
  const auto table_boundary = replace(input, "source = \"function\"\nexpression = \"sin(t)\"",
      "source = \"table\"\npath = \"boundary.csv\"\nt_column = \"t\"\nvalue_column = \"u\"");
  check(ConfigProcessor::parse(toml::parse(table_boundary)).boundaries.right.t_column == "t",
        "boundary table");
  const auto neumann = replace(input, "type = \"newton\"\nalpha = 2\nbeta = 0.5",
      "type = \"neumann\"");
  check(ConfigProcessor::parse(toml::parse(neumann)).boundaries.right.type == "neumann",
        "neumann boundary");

  rejects(replace(input, "config_version = 1", "config_version = 2"), "config_version");
  rejects(replace(input, "config_version = 1", "config_version = 1.0"), "invalid value type");
  rejects(replace(input, "config_version = 1", "config_version = 1\ntypo = 0"), "unknown");
  rejects(replace(input, "c = 1", "c = 0"), "c: must be > 0");
  rejects(replace(input, "c = 1", "c = nan"), "must be finite");
  rejects(replace(input, "L = 10.0", "sin = 10.0"), "reserved name");
  rejects(replace(input, "L = 10.0", "X = 10.0"), "reserved name");
  rejects(replace(input, "L = 10.0", "L = 10.0\nl = 2.0"), "differs only by case");
  rejects(replace(input, "L = 10.0", "\"1bad\" = 10.0"), "invalid or reserved");
  rejects(replace(input, "x = [0, 10]", "x = [10, 0]"), "upper bound");
  rejects(replace(input, "dx = 0.1", "dx = 0.3"), "multiple of step");
  rejects(replace(input, "dt = 0.05", "dt = 0"), "must be > 0");
  rejects(replace(input, "dt = 0.05", "dt = 0.2"), "CFL");
  const auto unstable = replace(replace(input, "dt = 0.05", "dt = 0.2"),
      "[solver]", "[solver]\ncheck_stability = false");
  check(!ConfigProcessor::parse(toml::parse(unstable)).solver.check_stability, "CFL switch");
  rejects(replace(unstable, "dt = 0.2", "dt = -0.2"), "must be > 0");
  rejects(replace(input, "[solver]", "[solver]\ncheck_stability = 1"), "invalid value type");
  rejects(replace(input, "explicit_central", "unknown"), "unsupported value");
  rejects(replace(input, "beta = 0.5", "beta = 0"), "beta: must be nonzero");
  rejects(replace(input, "beta = 0.5", ""), "missing required field");
  rejects(replace(input, "type = \"newton\"", "type = \"neumann\""), "inapplicable field");
  rejects(replace(input, "source = \"constant\"", "source = \"unknown\""), "unsupported value");
  rejects(replace(input, "value = 0", "value = true"), "expected number");
  rejects(replace(input, "value = 0", "value = 0\nexpression = \"x\""), "inapplicable field");
  rejects(replace(input, "expression = \"sin(t)\"", "expression = \"  \""), "must not be empty");
  rejects(replace(input, "expression = \"sin(t)\"", "expression = \"sin(unknown)\""),
          "unknown name");
  rejects(replace(input, "expression = \"sin(t)\"", "expression = \"max(t, 1)\""),
          "unknown name or function");
  rejects(replace(input, "expression = \"sin(t)\"", "expression = \"sin(\""),
          "expression");
  rejects(replace(input, "expression = \"sin(t)\"", "expression = \"x := 1\""),
          "unsupported character");
  rejects(replace(input, "expression = \"sin(t)\"", "expression = \"x^2\""),
          "unsupported character");
  rejects(replace(input, "expression = \"sin(t)\"", "expression = \"2x\""),
          "missing operator");
  rejects(replace(input, "expression = \"x - 5\"", "expression = \"sqrt(unknown)\""),
          "unknown name");
  const auto bad_domain = replace(input, "expression = \"sin(t)\"",
                                  "expression = \"sqrt(t - 1)\"");
  try {
    ConfigProcessor::parse(toml::parse(bad_domain)).boundaries.right.evaluate(0, 0);
    throw std::runtime_error("expected non-finite expression rejection");
  } catch (const std::runtime_error &error) {
    check(std::string(error.what()).find("finite") != std::string::npos,
          "non-finite expression diagnostic");
  }
  rejects(replace(input, "interval = [5, 10]", "interval = [6, 10]"), "gap");
  rejects(replace(input, "interval = [5, 10]", "interval = [4, 10]"), "overlap");
  rejects(replace(input, "interval = [5, 10]", "interval = [5, 9]"), "cover");
  rejects(replace(input, "expression = \"x - 5\"", "expression = \"x - 5\"\nvalue = 0"),
          "exactly one");
  rejects(replace(input, "x_column = \"x\"", "t_column = \"t\""), "x_column");
  rejects(replace(input, "[output]", "[output]\nsave_every = 0"), "save_every");
  rejects(replace(input, "[output]", "[output]\nfields = []"), "nonempty array");
  rejects(replace(input, "[output]", "[output]\nfields = [\"position\", \"position\"]"), "duplicate");
  rejects(replace(input, "[output]", "[output]\nfields = [\"energy\"]"), "unsupported value");

  // Проверка parse_file: пути берутся от файла, а не от рабочей директории.
  const auto file = std::filesystem::temp_directory_path() / "scienceutils-config-parser-test.toml";
  {
    std::ofstream stream(file);
    stream << input;
    check(static_cast<bool>(stream), "write test config");
  }
  const auto from_file = ConfigProcessor::parse_file(file);
  std::filesystem::remove(file);
  check(from_file.initials.velocity.path == (file.parent_path() / "data/velocity.csv").string(),
        "parse_file relative path");
  const auto sample = ConfigProcessor::parse_file("test.toml");
  check(std::abs(sample.initials.position.evaluate(10, 0) - 1.0) < 1e-12,
        "project sample expression");
  std::cout << "TOML parser tests passed\n";
}
