#pragma once

#include <array>
#include <cstddef>
#include <span>

template<typename T, int N>
requires (N > 0)
class ArraySpan {
private:
  std::span<std::array<T, N>, std::dynamic_extent> data;

public:
   template<size_t Extent>
   constexpr ArraySpan(std::span<std::array<T, N>, Extent> data)
     : data(data) {}

   constexpr ArraySpan(std::array<T, N>* data, size_t numVectors)
     : data(data, numVectors) {}

   constexpr std::array<T, N> operator[](int index) noexcept const {
     return coords[index];
   }

   T& operator[](int index) noexcept {
     return coords[index];
   }

   constexpr T operator[](int vectorIndex, int coordinateIndex) noexcept const {
     return data[vectorIndex][coordinateIndex];
   }

   T& operator[](int vectorIndex, int coordinateIndex) noexcept {
     return data[vectorIndex][coordinateIndex];
   }

   constexpr size_t size() const noexcept {
     return data.size();
   }

   static bool operator==(const ArraySpan<T, N>& a, const ArraySpan<T, N>& b) {
     if (a.size() != b.size()) {
       return false;
     }
     for (int i{ 0 }; i < a.size(); ++i) {
       for (int j{ 0 }; j < N; ++j) {
         if (a[i, j] != b[i, j]) {
           return false;
         }
       }
     }

     return true;
   }

};

namespace std {
  template <typename T, int N>
  struct hash<ArraySpan<T, N>> {
    size_t operator()(const ArraySpan<T, N>& span) const {
      size_t hash = 0;
      std::hash<T> hasher;
      for (int i{ 0 }; i < span.size(); ++i) {
        // TODO: Make this loop more "constexpr" using templates.
        for (int j{ 0 }; j < N; ++j) {
          hash ^= hasher(span[i, j]) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
      }
      return seed;
    }
  };
}
