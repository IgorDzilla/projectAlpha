#pragma once

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

/*
 * Абстрактный класс для любого Mesh-объекта.
 * */
template <std::size_t Dim> class Mesh {
  static_assert(Dim > 0, "Mesh dimension must be positive");

public:
  using Point = std::array<double, Dim>;
  using CellId = std::size_t;
  using NodeId = std::size_t;

  virtual ~Mesh() = default;

  // Общее количество ячеек
  virtual std::size_t cellCount() const noexcept = 0;

  // Координаты центра ячейки
  virtual Point cellCenter(CellId id) const = 0;

  // Длина, площадь или объем ячейки
  virtual double cellVolume(CellId id) const = 0;

  // Узлы — вершины ячеек, включая границу области. В 1D их N+1 для N ячеек.
  virtual std::size_t nodeCount() const noexcept = 0;
  virtual Point nodePosition(NodeId id) const = 0;
};

/*
 * Mesh в декартовой системе координат. Dim означает размерность.
 * */
template <std::size_t Dim> class CartesianMesh final : public Mesh<Dim> {

public:
  using Point = typename Mesh<Dim>::Point;
  using Index = std::array<std::size_t, Dim>;

private:
  Point origin_; // начало координат
  Point step_;   // шаг по каждой оси
  Point upperBound_;
  Index shape_;  // число ячеек на каждую ось
  Index strides_;
  Index nodeShape_;   // число узлов по каждой оси: shape[d] + 1
  Index nodeStrides_; // первая ось, как и у ячеек, меняется быстрее остальных

  std::size_t totalCells_ = 1;
  std::size_t totalNodes_ = 1;
  double volume_ = 1.0;

public:
  CartesianMesh(const Point &origin, const Point &step, const Index &shape)
      : origin_(origin), step_(step), shape_(shape) {
    for (std::size_t d = 0; d < Dim; ++d) {
      if (!std::isfinite(origin_[d]) || !std::isfinite(step_[d]) ||
          step_[d] <= 0.0 || shape_[d] == 0)
        throw std::invalid_argument("Mesh requires finite coordinates, positive steps and nonzero sizes");

      if (totalCells_ > std::numeric_limits<std::size_t>::max() / shape_[d])
        throw std::overflow_error("Mesh cell count is too large");

      if (shape_[d] == std::numeric_limits<std::size_t>::max())
        throw std::overflow_error("Mesh node count is too large");
      nodeShape_[d] = shape_[d] + 1;
      if (totalNodes_ > std::numeric_limits<std::size_t>::max() / nodeShape_[d])
        throw std::overflow_error("Mesh node count is too large");
      nodeStrides_[d] = totalNodes_;
      totalNodes_ *= nodeShape_[d];

      upperBound_[d] = std::fma(static_cast<double>(shape_[d]), step_[d], origin_[d]);
      if (!std::isfinite(upperBound_[d]) || upperBound_[d] <= origin_[d])
        throw std::invalid_argument("Mesh bounds cannot be represented");

      strides_[d] = totalCells_;
      totalCells_ *= shape_[d];
      volume_ *= step_[d];
      if (!std::isfinite(volume_) || volume_ <= 0.0)
        throw std::invalid_argument("Mesh cell volume cannot be represented");
    }
  }

  std::size_t cellCount() const noexcept override { return totalCells_; }

  std::size_t nodeCount() const noexcept override { return totalNodes_; }

  // В узле нет сдвига на половину шага, который используется для центра ячейки.
  Point nodePosition(std::size_t id) const override {
    const Index index = unflattenNode(id);
    Point result{};
    for (std::size_t d = 0; d < Dim; ++d)
      result[d] = std::fma(static_cast<double>(index[d]), step_[d], origin_[d]);
    return result;
  }

  /*
   * Вспомогательный метод. Возвращает центр ячейки.
   * */
  Point cellCenter(std::size_t id) const override {
    const Index index = unflatten(id);
    Point result{};

    for (std::size_t d = 0; d < Dim; ++d) {
      result[d] = std::fma(static_cast<double>(index[d]) + 0.5, step_[d], origin_[d]);
    }

    return result;
  }

  double cellVolume(std::size_t id) const override {
    if (id >= totalCells_)
      throw std::out_of_range("Cell id");
    return volume_;
  }

  // Координата грани по выбранной оси. Допустимы номера от 0 до shape[axis].
  double faceCoordinate(std::size_t axis, std::size_t faceIndex) const {
    if (axis >= Dim || faceIndex > shape_[axis])
      throw std::out_of_range("Mesh face index");
    return std::fma(static_cast<double>(faceIndex), step_[axis], origin_[axis]);
  }

  /*
   * Перевод многомерного индекса ЯЧЕЙКИ в её линейный номер.
   * */
  std::size_t flatten(const Index &index) const {
    std::size_t id = 0;

    for (std::size_t d = 0; d < Dim; ++d) {
      if (index[d] >= shape_[d])
        throw std::out_of_range("Mesh index");

      id += index[d] * strides_[d];
    }

    return id;
  }

  /*
   * Перевод линейного номера ЯЧЕЙКИ в многомерный индекс.
   * */
  Index unflatten(std::size_t id) const {
    if (id >= totalCells_)
      throw std::out_of_range("Cell id");

    Index index{};

    for (std::size_t d = 0; d < Dim; ++d) {
      index[d] = (id / strides_[d]) % shape_[d];
    }

    return index;
  }

  const Index &shape() const noexcept { return shape_; }

  // Индексы УЗЛОВ по каждой оси допустимы от 0 до shape[d] включительно.
  std::size_t flattenNode(const Index &index) const {
    std::size_t id = 0;
    for (std::size_t d = 0; d < Dim; ++d) {
      if (index[d] >= nodeShape_[d])
        throw std::out_of_range("Mesh node index");
      id += index[d] * nodeStrides_[d];
    }
    return id;
  }

  Index unflattenNode(std::size_t id) const {
    if (id >= totalNodes_)
      throw std::out_of_range("Node id");
    Index index{};
    for (std::size_t d = 0; d < Dim; ++d)
      index[d] = (id / nodeStrides_[d]) % nodeShape_[d];
    return index;
  }

  const Index &nodeShape() const noexcept { return nodeShape_; }

  const Point &step() const noexcept { return step_; }

  const Point &origin() const noexcept { return origin_; }

  const Point &upperBound() const noexcept { return upperBound_; }
};
