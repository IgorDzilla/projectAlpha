#include "toml-processor.hpp"

#include "exprtk.hpp"

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <mutex>
#include <stdexcept>
#include <string_view>

namespace tomlprocessor {
namespace {

bool identifier_start(char c) {
  return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || c == '_';
}

bool identifier_part(char c) {
  return identifier_start(c) || (c >= '0' && c <= '9');
}

// ExprTk поддерживает также присваивания, циклы и другие конструкции.
// Здесь пропускаются только токены языка формул из спецификации проекта.
std::string checked_formula(const std::string &input,
                            const std::map<std::string, double> &parameters) {
  enum class Last { start, value, operation, open, comma, function };
  Last last = Last::start;
  std::string result;
  for (std::size_t i = 0; i < input.size();) {
    const char c = input[i];
    if (std::isspace(static_cast<unsigned char>(c))) {
      ++i;
      continue;
    }
    if (std::isdigit(static_cast<unsigned char>(c)) || c == '.') {
      if (last == Last::value)
        throw std::runtime_error("missing operator before number");
      char *end = nullptr;
      const double value = std::strtod(input.c_str() + i, &end);
      if (end == input.c_str() + i || !std::isfinite(value))
        throw std::runtime_error("invalid or non-finite number");
      const auto length = static_cast<std::size_t>(end - (input.c_str() + i));
      result.append(input, i, length);
      i += length;
      last = Last::value;
    } else if (identifier_start(c)) {
      if (last == Last::value)
        throw std::runtime_error("missing operator before identifier");
      const auto begin = i++;
      while (i < input.size() && identifier_part(input[i]))
        ++i;
      const std::string name = input.substr(begin, i - begin);
      const bool function = name == "sin" || name == "cos" || name == "tan" ||
          name == "exp" || name == "log" || name == "sqrt" || name == "abs";
      if (function) {
        auto next = i;
        while (next < input.size() && std::isspace(static_cast<unsigned char>(input[next])))
          ++next;
        if (next == input.size() || input[next] != '(')
          throw std::runtime_error(name + ": expected function call");
        last = Last::function;
      } else if (name == "x" || name == "t" || name == "pi" || name == "e" ||
                 parameters.count(name)) {
        last = Last::value;
      } else {
        throw std::runtime_error(name + ": unknown name or function");
      }
      result += name;
    } else if (c == '(') {
      if (last == Last::value)
        throw std::runtime_error("missing operator before '('");
      result += c;
      last = Last::open;
      ++i;
    } else if (c == ')') {
      if (last != Last::value)
        throw std::runtime_error("unexpected ')'");
      result += c;
      last = Last::value;
      ++i;
    } else if (c == ',') {
      if (last != Last::value)
        throw std::runtime_error("unexpected ','");
      result += c;
      last = Last::comma;
      ++i;
    } else if (c == '+' || c == '-' || c == '*' || c == '/') {
      if ((c == '*' || c == '/') && last != Last::value)
        throw std::runtime_error("unexpected arithmetic operator");
      if (c == '*' && i + 1 < input.size() && input[i + 1] == '*') {
        // В конфигурации степень записывается как **, в ExprTk — как ^.
        result += '^';
        i += 2;
      } else {
        result += c;
        ++i;
      }
      last = Last::operation;
    } else {
      throw std::runtime_error(std::string("unsupported character '") + c + "'");
    }
  }
  if (last != Last::value)
    throw std::runtime_error("incomplete expression");
  return result;
}

} // namespace

struct MathExpression::Impl {
  double x = 0.0;
  double t = 0.0;
  std::map<std::string, double> parameters;
  exprtk::symbol_table<double> symbols;
  exprtk::expression<double> expression;
  std::mutex mutex;

  Impl(const std::string &formula, const std::map<std::string, double> &values)
      : parameters(values) {
    // ExprTk хранит ссылки на переменные, поэтому их адреса должны жить
    // столько же, сколько скомпилированное выражение.
    if (!symbols.add_variable("x", x) || !symbols.add_variable("t", t) ||
        !symbols.add_constant("pi", std::acos(-1.0)) ||
        !symbols.add_constant("e", std::exp(1.0)))
      throw std::runtime_error("cannot register expression coordinates");
    for (auto &[name, value] : parameters)
      if (!symbols.add_constant(name, value))
        throw std::runtime_error(name + ": cannot register parameter");
    expression.register_symbol_table(symbols);

    exprtk::parser<double> parser;
    using Settings = exprtk::parser<double>::settings_store;
    parser.settings().disable_all_base_functions()
        .enable_base_function(Settings::e_bf_sin)
        .enable_base_function(Settings::e_bf_cos)
        .enable_base_function(Settings::e_bf_tan)
        .enable_base_function(Settings::e_bf_exp)
        .enable_base_function(Settings::e_bf_log)
        .enable_base_function(Settings::e_bf_sqrt)
        .enable_base_function(Settings::e_bf_abs)
        .disable_all_control_structures()
        .disable_all_logic_ops()
        .disable_all_assignment_ops()
        .disable_all_inequality_ops()
        .disable_local_vardef();
    if (!parser.compile(checked_formula(formula, parameters), expression))
      throw std::runtime_error("invalid expression: " + parser.error());
  }
};

MathExpression::MathExpression(const std::string &formula,
                               const std::map<std::string, double> &parameters)
    : impl_(std::make_shared<Impl>(formula, parameters)) {}

double MathExpression::evaluate(double x, double t) const {
  if (!std::isfinite(x) || !std::isfinite(t))
    throw std::runtime_error("expression coordinates must be finite");
  std::lock_guard<std::mutex> lock(impl_->mutex);
  impl_->x = x;
  impl_->t = t;
  const double value = impl_->expression.value();
  if (!std::isfinite(value))
    throw std::runtime_error("expression result must be finite");
  return value;
}

double Segment::evaluate(double x, double t) const {
  if (value)
    return *value;
  if (compiled_expression)
    return compiled_expression->evaluate(x, t);
  throw std::runtime_error("segment has no compiled expression");
}

double SourceFunction::evaluate(double x, double t) const {
  if (source == "constant")
    return value;
  if (source == "function" && compiled_expression)
    return compiled_expression->evaluate(x, t);
  throw std::runtime_error("source '" + source + "' cannot be evaluated directly");
}

} // namespace tomlprocessor
