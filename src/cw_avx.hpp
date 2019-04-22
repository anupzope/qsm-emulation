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

namespace qsm {

/**
 * Specialization of bsm::cw_traits<N> for a control word type of 3 bits.
 */
template<>
struct cw_traits<3> {
  using type = uint8_t;
  static constexpr int NParts = 8;
};

/**
 * Specialization of bsm::cw_traits<N> for a control word type of 4 bits.
 */
template<>
struct cw_traits<4> {
  using type = uint16_t;
  static constexpr int NParts = 16;
};

/**
 * Specialization of bsm::cw_traits<N> for a control word type of 5 bits.
 */
template<>
struct cw_traits<5> {
  using type = uint32_t;
  static constexpr int NParts = 32;
};

/**
 * Specialization of bsm::cw_traits<N> for a control word type of 6 bits.
 */
template<>
struct cw_traits<6> {
  using type = uint64_t;
  static constexpr int NParts = 64;
};

/**
 * Specialization of bsm::cw_traits<N> for a control word type of 7 bits.
 */
template<>
struct cw_traits<7> {
  using type = __m128i;
  static constexpr int NParts = 128;
};

/**
 * Specialization of bsm::cw_traits<N> for a control word type of 8 bits.
 */
template<>
struct cw_traits<8> {
  using type = __m256i;
  static constexpr int NParts = 256;
};

/**
 * Specialization of bsm::cwiterator<N> for a control word of 3 bits.
 */
class cwiterator<3> {
public:
  cwiterator(typename cw_traits<3>::type v) : value(v) {}
  
  bool has_next() const {
    return testnotzero(value);
  }
  
  int next() {
    int index = ls1bindex(value);
    value &= (value - 1);
    return index;
  }
  
private:
  typename cw_traits<3>::type value;
};

/**
 * Specialization of bsm::cwiterator<N> for a control word of 4 bits.
 */
class cwiterator<4> {
public:
  cwiterator(typename cw_traits<4>::type v) : value(v) {}
  
  bool has_next() const {
    return testnotzero(value);
  }
  
  int next() {
    int index = ls1bindex(value);
    value &= (value - 1);
    return index;
  }
  
private:
  typename cw_traits<4>::type value;
};

/**
 * Specialization of bsm::cwiterator<N> for a control word of 5 bits.
 */
class cwiterator<5> {
public:
  cwiterator(typename cw_traits<5>::type v) : value(v) {}
  
  bool has_next() const {
    return testnotzero(value);
  }
  
  int next() {
    int index = ls1bindex(value);
    value &= (value - 1);
    return index;
  }
  
private:
  typename cw_traits<5>::type value;
};

/**
 * Specialization of bsm::cwiterator<N> for a control word of 6 bits.
 */
class cwiterator<6> {
public:
  cwiterator(typename cw_traits<6>::type v) : value(v) {}
  
  bool has_next() const {
    return testnotzero(value);
  }
  
  int next() {
    int index = ls1bindex(value);
    value &= (value - 1);
    return index;
  }
  
private:
  typename cw_traits<6>::type value;
};

/**
 * Specialization of bsm::cwiterator<N> for a control word of 7 bits.
 */
template<>
class cwiterator<7> {
public:
  cwiterator(typename cw_traits<7>::type v) : value(v), idx(0), offset(0) {
  }
  
  /** Tests if there is any set bit in the cw associated with this iterator. */
  bool has_next() const {
    return testnotzero(value);
  }
  
  /** Returns the index of next set bit in the cw associated with this iterator. */
  int next() {
    while(idx < 4 && !values[idx]) {
      idx++; offset += 32;
    }
    int index = ls1bindex(values[idx]) + offset;
    values[idx] &= (values[idx] - 1);
    return index;
  }
  
private:
  int offset;
  int idx;
  union {
    uint32_t values[4];
    typename cw_traits<7>::type value;
  };
};

/**
 * Specialization of bsm::cwiterator<N> for a control word of 8 bits.
 */
template<>
class cwiterator<8> {
public:
  cwiterator(typename cw_traits<8>::type v) : value(v), idx(0), offset(0) {
  }
  
  /** Tests if there is any set bit in the cw associated with this iterator. */
  bool has_next() const {
    return testnotzero(value);
  }
  
  /** Returns the index of next set bit in the cw associated with this iterator. */
  int next() {
    while(idx < 8 && !values[idx]) {
      idx++; offset += 32;
    }
    int index = ls1bindex(values[idx]) + offset;
    values[idx] &= (values[idx] - 1);
    return index;
  }
  
private:
  int offset;
  int idx;
  union {
    uint32_t values[8];
    typename cw_traits<8>::type value;
  };
};

} // end: namespace qsm
