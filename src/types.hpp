//##############################################################################
//# Copyright 2019, Mississippi State University
//# 
//# This file is part of qsm-emulation.
//# 
//# qsm-emulation is free software: you can redistribute it and/or modify
//# it under the terms of the GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//# 
//# qsm-emulation is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//# 
//# You should have received a copy of the GNU General Public License
//# along with this program.  If not, see <https://www.gnu.org/licenses/>.
//##############################################################################

#pragma once

#include <array>

template<typename T>
inline T reciprocal(T const & v) {
  return (T)1/v;
}

template<typename T, int R, int C>
class matrix;

template<typename T, int R, int C>
inline std::ostream & operator<<(std::ostream & stream, matrix<T, R, C> const & val) {
  for(int r = 0; r < R; ++r) {
    for(int c = 0; c < C; ++c) {
      stream << val(r, c) << " ";
    }
  }
  return stream;
}

template<typename T, int R, int CR, int C>
inline matrix<T, R, C> operator*(matrix<T, R, CR> const & lhs, matrix<T, CR, C> const & rhs) {
  matrix<T, R, C> res;
  for(int r = 0; r < R; ++r) {
    for(int c = 0; c < C; ++c) {
      T sum = 0;
      for(int k = 0; k < CR; ++k) {
        sum += lhs(r, k) * rhs(k, c);
      }
      res(r, c) = sum;
    }
  }
  
  return res;
}

template<typename T, int R, int C>
class matrix {
  T e[R][C];
  
public:
  matrix() {
  }
  
  matrix(T val) {
    for(int r = 0; r < R; ++r) {
      for(int c = 0; c < C; ++c) {
        e[r][c] = val;
      }
    }
  }
  
  matrix<T, R, C> & operator+=(matrix<T, R, C> const & rhs) {
    for(int r = 0; r < R; ++r) {
      for(int c = 0; c < C; ++c) {
        e[r][c] += rhs.e[r][c];
      }
    }
    return *this;
  }
  
  matrix<T, R, C> & operator-=(matrix<T, R, C> const & rhs) {
    for(int r = 0; r < R; ++r) {
      for(int c = 0; c < C; ++c) {
        e[r][c] -= rhs.e[r][c];
      }
    }
    return *this;
  }
  
  T & operator()(int r, int c) {
    return e[r][c];
  }
  
  T const & operator()(int r, int c) const {
    return e[r][c];
  }
};

template<typename T, int N>
class vec_t {
private:
  T m_data[N];
  
public:
  vec_t() {
  }
  
  vec_t(T const & scalar) {
    for(int i = 0; i < N; ++i)
      m_data[i] = scalar;
  }
  
  vec_t(vec_t<T, N> const & src) {
    for(int i = 0; i < N; ++i)
      m_data[i] = src.m_data[i];
  }
  
  template<typename... V>
  vec_t(V... args) : m_data{args...} {
  }
  
  vec_t<T, N> & operator=(T const & scalar) {
    for(int i = 0; i < N; ++i)
      m_data[i] = scalar;
    return *this;
  }
  
  vec_t<T, N> & operator=(vec_t<T, N> const & rhs) {
    for(int i = 0; i < N; ++i)
      m_data[i] = rhs.m_data[i];
    return *this;
  }
  
  vec_t<T, N> & operator=(std::initializer_list<T> l) {
    auto iter = l.begin();
    for(int i = 0; i < N; ++i)
      m_data[i] = *iter++;
    return *this;
  }
  
  vec_t<T, N> & operator+=(vec_t<T, N> const & rhs) {
    for(int i = 0; i < N; ++i)
      m_data[i] += rhs.m_data[i];
    return *this;
  }
  
  vec_t<T, N> & operator-=(vec_t<T, N> const & rhs) {
    for(int i = 0; i < N; ++i)
      m_data[i] -= rhs.m_data[i];
    return *this;
  }
  
  vec_t<T, N> & operator*=(T const & scalar) {
    for(int i = 0; i < N; ++i)
      m_data[i] *= scalar;
    return *this;
  }
  
  vec_t<T, N> & operator/=(T const & scalar) {
    T const rscalar  = reciprocal(scalar);
    for(int i = 0; i < N; ++i)
      m_data[i] *= rscalar;
    return *this;
  }
  
  T & operator[](int const index) {
    return m_data[index];
  }
  
  T const & operator[](int const index) const {
    return m_data[index];
  }
};

template<typename T, int N>
inline vec_t<T, N> operator+(vec_t<T, N> const & lhs, vec_t<T, N> const & rhs) {
  vec_t<T, N> res;
  for(int i = 0; i < N; ++i)
    res[i] = lhs[i] + rhs[i];
  return res;
}

template<typename T, int N>
inline vec_t<T, N> operator-(vec_t<T, N> const & lhs, vec_t<T, N> const & rhs) {
  vec_t<T, N> res;
  for(int i = 0; i < N; ++i)
    res[i] = lhs[i] - rhs[i];
  return res;
}

template<typename T, int N>
inline vec_t<T, N> subtract_reverse(vec_t<T, N> const & lhs, vec_t<T, N> const & rhs) {
  vec_t<T, N> res;
  for(int i = N-1; i >= 0; --i)
    res[i] = lhs[i] - rhs[i];
  return res;
}

template<typename T, int N>
inline vec_t<T, N> operator*(vec_t<T, N> const & lhs, T const & scalar) {
  vec_t<T, N> res;
  for(int i = 0; i < N; ++i)
    res[i] = lhs[i]*scalar;
  return res;
}

template<typename T, int N>
inline vec_t<T, N> operator*(T const & scalar, vec_t<T, N> const & rhs) {
  vec_t<T, N> res;
  for(int i = 0; i < N; ++i)
    res[i] = rhs[i]*scalar;
  return res;
}

template<typename T, int N>
inline vec_t<T, N> operator/(vec_t<T, N> const & lhs, T const & scalar) {
  T const rscalar  = reciprocal(scalar);
  vec_t<T, N> res;
  for(int i = 0; i < N; ++i)
    res[i] = lhs[i]*rscalar;
  return res;
}

template<typename T>
inline vec_t<T, 3> cross(vec_t<T, 3> const & lhs, vec_t<T, 3> const & rhs) {
  //vec_t<T, 3> res;
  //res[0] = lhs[1]*rhs[2]-lhs[2]*rhs[1];
  //res[1] = lhs[2]*rhs[0]-lhs[0]*rhs[2];
  //res[2] = lhs[0]*rhs[1]-lhs[1]*rhs[0];
  //return res;
  return {
    lhs[1]*rhs[2]-lhs[2]*rhs[1],
    lhs[2]*rhs[0]-lhs[0]*rhs[2],
    lhs[0]*rhs[1]-lhs[1]*rhs[0]
  };
}

template<typename T, int N>
inline T dot(vec_t<T, N> const & lhs, vec_t<T, N> const & rhs) {
  T res = (T)0;
  for(int i = 0; i < N; ++i)
    res += lhs[i]*rhs[i];
  return res;
}

template<typename T, int N>
inline T dot_reverse(vec_t<T, N> const & lhs, vec_t<T, N> const & rhs) {
  T res = (T)0;
  for(int i = N-1; i >= 0; --i)
    res += lhs[i]*rhs[i];
  return res;
}

template<typename ... T>
class pack_store {
  std::tuple<T * ...> ptrs;
  
  template<size_t ... I>
  void allocate1(size_t const size, std::index_sequence<I...>) {
    ptrs = { new T[size] ... };
  }
  
  template<size_t ... I>
  void deallocate1(std::index_sequence<I...>) {
    using swallow = int[];
    (void)swallow{0, (delete[] std::get<I>(ptrs), 0)... };
  }
  
public:
  void allocate(size_t const size) {
    allocate1(size, std::make_index_sequence<sizeof...(T)>{});
  }
  
  void deallocate() {
    deallocate1(std::make_index_sequence<sizeof...(T)>{});
  }
  
  std::tuple<T * ...> get() {
    return ptrs;
  }
  
  std::tuple<T const * ...> get_const() const {
    return std::tuple<T const * ...>(ptrs);
  }
};
