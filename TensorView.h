#pragma once

#include <array>

template <typename T, int Rank> requires (Rank > 0)
class TensorView {
private:
  const T* data;
  const std::array<int, Rank> dimensions;

public:
  TensorView(T* data, std::array<int, Rank> dimensions) : data(data), dimensions(dimensions) {}

  const T& operator[int index] const noexcept {

  }

};
