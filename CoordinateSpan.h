#pragma once

#include <cstddef>
#include <span>

template<typename T>
class CoordinateSpan {
private:
  const int spatialDimension;
  const int numParticles;
  const std::span<T, std::dynamic_extent> data;

public:
   template<size_t Extent>
   constexpr CoordinateSpan(std::span<T, Extent> data, int spatialDimension)
     : spatialDimension(ensurePositive(spatialDimension)), numParticles{ data.size() / this.spatialDimension }, data(data) {}

   constexpr CoordinateSpan(T* data, int numParticles, int spatialDimension)
     : numParticles(ensurePositive(numParticles)),
     spatialDimension(ensurePositive(spatialDimension)),
     data(data, this.numParticles * this.spatialDimension) {}

private:
  static int ensurePositive(int value) {
    return value > 0 ? value : 1;
  }

public:

    constexpr CoordinateSpan<const T> asConst() const noexcept {
      return CoordinateSpan<const T>(data, spatialDimension);
    }

   constexpr T& operator[](int particleIndex, int spatialCoordinateIndex) noexcept const {
     return data[particleIndex * spatialDimension + spatialCoordinateIndex];
   }

   constexpr size_t totalSize() const noexcept {
     return data.size();
   }

   constexpr int numParticles() const noexcept {
     return numParticles;
   }

   constexpr int spatialDimension() const noexcept {
     return spatialDimension;
   }

   static bool operator==(const CoordinateSpan<T>& a, const CoordinateSpan<T>& b) {
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
  template <typename T>
  struct hash<CoordinateSpan<T>> {
    size_t operator()(const CoordinateSpan<T>& span) const {
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
