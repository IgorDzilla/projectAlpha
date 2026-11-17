#include "EulerConfig.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <set>
#include <stdexcept>
#include <string_view>

namespace {
[[noreturn]] void fail(const std::string& path, const std::string& message) {
  throw std::invalid_argument(path + ": " + message);
}
void keys(const toml::table& table, const std::string& path,
          std::initializer_list<std::string_view> allowed) {
  for (const auto& [key, node] : table) {
    (void)node;
    if (std::find(allowed.begin(), allowed.end(), key.str()) == allowed.end())
      fail(path + "." + std::string(key.str()), "unknown or inapplicable field");
  }
}
const toml::node& required(const toml::table& t, std::string_view key, const std::string& path) {
  const auto* node = t.get(key);
  if (!node) fail(path + "." + std::string(key), "required field missing");
  return *node;
}
const toml::table& table(const toml::table& t, std::string_view key, const std::string& path) {
  const auto* value = required(t, key, path).as_table();
  if (!value) fail(path + "." + std::string(key), "expected table");
  return *value;
}
std::string string(const toml::node& n, const std::string& path) {
  auto value = n.value_exact<std::string>();
  if (!value || value->find_first_not_of(" \t\r\n") == std::string::npos ||
      value->find('\0') != std::string::npos)
    fail(path, "expected nonempty string without NUL");
  return *value;
}
void choice(const toml::table& t, std::string_view key, const std::string& path,
            const std::string& expected) {
  const auto value = string(required(t, key, path), path + "." + std::string(key));
  if (value != expected)
    fail(path + "." + std::string(key), "unsupported '" + value + "'; expected '" + expected + "'");
}
double number(const toml::node& n, const std::string& path) {
  if (!n.is_number()) fail(path, "expected number");
  const double v = n.is_integer() ? static_cast<double>(*n.value_exact<std::int64_t>())
                                  : *n.value_exact<double>();
  if (!std::isfinite(v)) fail(path, "must be finite");
  return v;
}
double number(const toml::table& t, std::string_view key, const std::string& path) {
  return number(required(t, key, path), path + "." + std::string(key));
}
std::size_t positiveInteger(const toml::node& n, const std::string& path) {
  auto v = n.value_exact<std::int64_t>();
  if (!v || *v <= 0 || static_cast<std::uint64_t>(*v) >= std::numeric_limits<std::size_t>::max())
    fail(path, "expected positive representable integer");
  return static_cast<std::size_t>(*v);
}
std::pair<double,double> interval(const toml::table& t, std::string_view key, bool allowEqual) {
  const std::string path = "domain." + std::string(key);
  const auto* array = required(t, key, "domain").as_array();
  if (!array || array->size() != 2) fail(path, "expected two numbers");
  const double a = number(*array->get(0), path), b = number(*array->get(1), path);
  if (b < a || (!allowEqual && b == a) || !std::isfinite(b-a))
    fail(path, "invalid interval");
  return {a,b};
}
EulerPrimitiveState primitive(const toml::table& t, const std::string& path, const EulerEquations& gas) {
  keys(t, path, {"density", "velocity", "pressure"});
  EulerPrimitiveState w{number(t,"density",path), number(t,"velocity",path), number(t,"pressure",path)};
  try { gas.toConservative(w); gas.soundSpeed(w); }
  catch (const std::exception& e) { fail(path,e.what()); }
  return w;
}
}

EulerConfig parseEulerConfig(const toml::table& root, const std::filesystem::path& baseDirectory) {
  keys(root,"root",{"config_version","problem","parameters","domain","grid","initials","solver","output"});
  const auto version = required(root,"config_version","root").value_exact<std::int64_t>();
  if (!version || *version != 3)
    fail("config_version", "this CLI requires exact Euler Riemann TOML v3");
  const auto& problem = table(root,"problem","root");
  keys(problem,"problem",{"type"}); choice(problem,"type","problem","euler_1d");
  EulerConfig c;
  const auto& parameters = table(root,"parameters","root");
  keys(parameters,"parameters",{"gamma"});
  c.gamma = number(parameters,"gamma","parameters");
  if (c.gamma <= 1.0) fail("parameters.gamma","must exceed one");
  const EulerEquations gas(c.gamma);
  const auto& domain = table(root,"domain","root");
  keys(domain,"domain",{"x","t"});
  const auto x = interval(domain,"x",false), t = interval(domain,"t",true);
  c.xMin=x.first; c.xMax=x.second; c.tStart=t.first; c.tEnd=t.second;
  const auto& grid = table(root,"grid","root");
  keys(grid,"grid",{"points"});
  c.points = positiveInteger(required(grid,"points","grid"),"grid.points");
  const double dx = (c.xMax-c.xMin)/static_cast<double>(c.points);
  if (!(dx > 0.0) || c.xMin + dx == c.xMin || c.xMax - dx == c.xMax)
    fail("grid.points","sampling spacing cannot be represented on this domain");
  const auto& initials = table(root,"initials","root");
  keys(initials,"initials",{"type","discontinuity","left","right"});
  choice(initials,"type","initials","riemann");
  c.discontinuity = number(initials,"discontinuity","initials");
  if (!(c.discontinuity > c.xMin && c.discontinuity < c.xMax))
    fail("initials.discontinuity","must lie strictly inside domain.x");
  c.left=primitive(table(initials,"left","initials"),"initials.left",gas);
  c.right=primitive(table(initials,"right","initials"),"initials.right",gas);
  const auto& solver=table(root,"solver","root");
  keys(solver,"solver",{"method","exact"});
  choice(solver,"method","solver","exact_riemann");
  if (solver.contains("exact")) {
    const auto& exact=table(solver,"exact","solver");
    keys(exact,"solver.exact",{"pressure_rtol","residual_rtol","max_iterations","max_bracket_expansions"});
    if (const auto* n=exact.get("pressure_rtol")) c.exact.pressureRtol=number(*n,"solver.exact.pressure_rtol");
    if (const auto* n=exact.get("residual_rtol")) c.exact.residualRtol=number(*n,"solver.exact.residual_rtol");
    if (const auto* n=exact.get("max_iterations")) c.exact.maxIterations=positiveInteger(*n,"solver.exact.max_iterations");
    if (const auto* n=exact.get("max_bracket_expansions")) c.exact.maxBracketExpansions=positiveInteger(*n,"solver.exact.max_bracket_expansions");
    if (c.exact.pressureRtol <= 0.0 || c.exact.residualRtol <= 0.0)
      fail("solver.exact","tolerances must be positive");
  }
  if (root.contains("output")) {
    const auto& out=table(root,"output","root");
    keys(out,"output",{"directory","format","fields"});
    choice(out,"format","output","csv");
    EulerOutputConfig output;
    const auto directory = string(required(out,"directory","output"),"output.directory");
    output.directory=std::filesystem::absolute(baseDirectory / directory).lexically_normal();
    if (out.contains("fields")) {
      const auto* fields=required(out,"fields","output").as_array();
      if (!fields || fields->empty()) fail("output.fields","expected nonempty array");
      const std::set<std::string> allowed{"density","momentum","energy","velocity","pressure"};
      std::set<std::string> seen;
      output.fields.clear();
      for (const auto& n : *fields) {
        const auto name=string(n,"output.fields");
        if (!allowed.count(name) || !seen.insert(name).second)
          fail("output.fields","unknown or duplicate field '"+name+"'");
        output.fields.push_back(name);
      }
    }
    c.output=std::move(output);
  }
  return c;
}

EulerConfig readEulerConfig(const std::filesystem::path& path) {
  const auto absolute=std::filesystem::absolute(path).lexically_normal();
  return parseEulerConfig(toml::parse_file(absolute.string()),absolute.parent_path());
}
