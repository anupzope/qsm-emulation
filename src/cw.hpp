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

#include "bitmanip.hpp"

namespace qsm {

/**
 * Traits for a control word that can be used to process maximum 
 * of \c N bits.
 * 
 * For each target platform, specializations of this class must 
 * be provided. A secialization of this class must provide 
 * following members.
 * -# A \c typedef that defines the data type used to store a control 
 *    word that processes maximum of \c N bits.
 *    \code
 *      typedef <c++type> type;
 *    \endcode
 * -# A \c static \c constexpr field of type \c int that specifies 
 *    the maximum number of partitions a control word of \c N bits 
 *    represents. This is generally set to \f$ 2^N \f$.
 *    \code
 *      static constexpr int NParts = <value-of-2^N>;
 *    \endcode
 *
 * \tparam N Maximum number of bits processed by a control word 
 *         represented by this trait.
 */
template<int N>
struct cw_traits {
};

/**
 * Iterator for a control word of \c N bits.
 * 
 * For each target platform, specializations of this class must 
 * be provided. The iterator for a control word may be 
 * destructive, that is, it may change the value of the control 
 * word. Therefore, the iterator always works on a copy of the 
 * control word. A specialization of this class must provide 
 * following members.
 * -# Ctor that takes a control word and makes a copy of it 
 *    for processing:
 *    \code
 *      cwiterator(typename cw_traits<N>::type v);
 *    \endcode
 * -# Method to check if there is any set bit left in the control 
 *    word.
 *    \code
 *      bool has_next() const;
 *    \endcode
 * -# Method to get the index of next set bit in the control word.
 *    \code
 *      int next();
 *    \endcode
 */
template<int N>
class cwiterator {
};

} // end: namespace qsm

// Include platform specific implementations of cw_traits and 
// cwiterator classes.
#if defined(__AVX__)
#include "cw_avx.hpp"
#else
#error "Platform specific implementation of bitmanip.hpp not found"
#endif

namespace qsm {

/**
 * Type for a control word of \c N bits.
 */
template<int N>
class cw {
public:
  /** Type fo iterator for this cw type. */
  using iterator = cwiterator<N>;
  
public:
  /** Default ctor. */
  cw() {}
  
  /** Ctor that takes a value of the control word. */
  cw(typename cw_traits<N>::type value) : value(value) {}
  
  /** Starts iteration over 'on' bits of this cw. */
  iterator begin() const {
    return iterator(value);
  }
  
  /** Tests if this cw has all the bits set to 'off'. */
  bool empty() const {
    return testzero(value);
  }
  
  /** Sets all the bits of this cw to 'off'. */
  void clear() {
    value = setzero<typename cw_traits<N>::type>();
  }
  
  bool is_set(int index) const {
    return value | (cw_traits<N>::type)((cw_traits<N>::type)0x1u << index);
  }
  
  /** Sets bit at index \c index in this cw to 'on'. */
  void set(int index) {
    value = setbit(value, index);
  }
  
  /** Tests if this control word and \c rhs are equal. */
  bool operator==(const cw<N> &rhs) const {
    return testequal(value, rhs.value);
  }
  
  /** Tests of this control word and \c rhs are not equal. */
  bool operator!=(const cw<N> &rhs) const {
    return testnotequal(value, rhs.value);
  }
  
  /** Returns the number of 'on' bits in this cw. */
  int count() const {
    return bit1count(value);
  }
  
private:
  typename cw_traits<N>::type value;
};

} // namespace qsm
