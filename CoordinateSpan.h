#pragma once

#include <assert.h>
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
     : spatialDimension(ensurePositive(spatialDimension)), numParticles{ data.size() / this->spatialDimension }, data(data) {}

   constexpr CoordinateSpan(T* data, int numParticles, int spatialDimension)
     : numParticles(ensurePositive(numParticles)),
     spatialDimension(ensurePositive(spatialDimension)),
     data(data, this->numParticles * this->spatialDimension) {}

private:
  static int ensurePositive(int value) {
    return value > 0 ? value : 1;
  }

public:

    constexpr CoordinateSpan<const T> asConst() const noexcept {
      return CoordinateSpan<const T>(data, spatialDimension);
    }

    constexpr T* raw() noexcept {
      return data.data();
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

   constexpr CoordinateSpan<T>& operator+=(const CoordinateSpan<T>& rhs) {
     assert(rhs.numParticles == numParticles);
     assert(rhs.spatialDimension == spatialDimension);
     
     int totalSize{ totalSize() };
     T *data = this->data.data();
     const T* rhsData = rhs.data.data();
     for (int i{ 0 }; i < totalSize; ++i) {
       data[i] += rhsData[i];
     }

     return *this;
   }

   constexpr CoordinateSpan<T>& operator-=(const CoordinateSpan<T>& rhs) {
     assert(rhs.numParticles == numParticles);
     assert(rhs.spatialDimension == spatialDimension);

     int totalSize{ totalSize() };
     T *data = this->data.data();
     const T* rhsData = rhs.data.data();
     for (int i{ 0 }; i < totalSize; ++i) {
       data[i] -= rhsData[i];
     }

     return *this;
   }

   constexpr CoordinateSpan<T>& operator*=(double rhs) {
     int totalSize{ totalSize() };
     T *data = this->data.data();
     for (int i{ 0 }; i < totalSize; ++i) {
       data[i] *= rhs;
     }

     return *this;
   }

   constexpr CoordinateSpan<T>& operator-() {
      int totalSize{ totalSize() };
     T *data = this->data.data();
     for (int i{ 0 }; i < totalSize; ++i) {
       data[i] = -data[i];
     }
     return *this;
   }

   static bool operator==(const CoordinateSpan<T>& a, const CoordinateSpan<T>& b) {
     if (a.size() != b.size()) {
       return false;
     }
     for (int i{ 0 }; i < a.size(); ++i) {
       for (int j{ 0 }; j < a.spatialDimension(); ++j) {
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
        for (int j{ 0 }; j < span.spatialDimension(); ++j) {
          hash ^= hasher(span[i, j]) + 0x9e3779b9 + (hash << 6) + (hash >> 2);
        }
      }
      return hash;
    }
  };
}
