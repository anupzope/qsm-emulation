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

#include <machineparams.hpp>
#include <utility.hpp>

namespace qsm {

template<typename W>
class cacheline_t;

template<typename W>
class mask_t;

/**
 * Abstraction of a mask for a cacheline.
 * 
 * The stream in and stream out operations of the class cacheline_t 
 * require masking to accurately control which words in the cacheline 
 * are written and read from the DRAM. This class provides the a way 
 * to represent the mask. Since different processors have different 
 * ways to perform masking on the cachelines. An implementation of 
 * this class is particular to a specific processor. Therefore, this 
 * class needs to be specialized for each target processor. Each specialization 
 * must provide implementations of the methods of this class.
 * 
 * \tparam W word type
 */
template<typename W>
class mask_t {
public:
  /** Default constructor. Creates a mask with all bits unset. */
  mask_t() {
    static_assert(false, "Not implemented");
  }
  
  /**
   * Creates a mask with unset bits at <tt>[0, step)</tt> and set bits at 
   * <tt>[step, end)</tt>. \c step is zero based.
   */
  mask_t(int step) {
    static_assert(false, "Not implemented");
  }
  
  /** Sets all bits of this mask. */
  void set() {
    static_assert(false, "Not implemented");
  }
  
  /** resets all bits in the mask. */
  void reset() {
    static_assert(false, "Not implemented");
  }
  
  /**
   * Unsets bits at [0, step) and sets bits at [step, end) of this mask. 
   * \c step is zero based.
   */
  void stepup(int step) {
    static_assert(false, "Not implemented");
  }
  
  /**
   * Sets bits at [0, step) and unsets bits at [step, end) of this mask. 
   * \c step is zero based.
   */
  void stepdown(int step) {
    static_assert(false, "Not implemented");
  }
  
  /** Flips the bits of this mask. */
  void flip() {
    static_assert(false, "Not implemented");
  }
  
  /** Tests if this mask is zero. */
  bool is_zero() const {
    static_assert(false, "Not implemented");
  }
  
  /** Xors this mask with rhs. */
  mask_t<W> & operator^=(mask_t<W> const & rhs) {
    static_assert(false, "Not implemented");
  }
};

/**
 * Abstraction of a cacheline.
 *
 * The data transfers between the DRAM and a core in chunks called cachelines. 
 * This class abstracts a cacheline. It consists of a sequence of a fixed 
 * number of words. On most of the hardware platforms, one byte is one word. 
 * But for software abstraction, we can have any C++ primitive type as a word 
 * e.g. char, short, int etc. The template parameter \c W specifies the type 
 * of word. This class needs to be specialized for each target processor. 
 * Often, the cachelines are aligned on certain address boundaries. The 
 * specializations must specify the alignment requirement using the \c alignas 
 * clause. In addition, they must provide implementations of the methods 
 * described in this skeleton. In addition to this, the specializations must 
 * also provide definition of following <tt>static constexpr</tt>s
 * -# <tt>int WordSize</tt> - the number of words per cacheline
 * -# <tt>uintptr_t ByteSize</tt> - number by which a <tt>void*</tt> converted 
 *    needs to advance into the next cache line
 * 
 * They also need to implement following static methods that round down and up 
 * a given pointer on cacheline boundaries.
 * <tt>template<typename T> static T * round_down(T * ptr)</tt>
 * <tt>template<typename T> static T * round_up(T * ptr)</tt>
 * 
 * \tparam W word type
 */
template<typename W>
class cacheline_t {
public:
  /**
   * Writes the words in this cacheline to the global memory in consecutive 
   * locations starting at an address \c memptr. This method must perform 
   * a direct write, bypassing the cache coherency protocol. \c memptr must be 
   * cacheline aligned.
   */
  void stream_out(W * memptr) {
    static_assert(false, "Not implemented");
  }
  
  void stream_out(cacheline_t<W> * memptr) {
    static_assert(false, "Not implemented");
  }
  
  /**
   * Writes the words in this cacheline masked by a mask to the global memory 
   * in consecutive locations starting at an address \c memptr. This method 
   * must perform a direct write, bypassing the cache coherency protocol. 
   * Only those words in this cacheline that are at positions specified as 'on' 
   * in the mask are written to global memory. \c memptr must be cacheline 
   * aligned.
   */
  void masked_stream_out(W * memptr, mask_t<W> & mask) {
    static_assert(false, "Not implemented");
  }
  
  void masked_stream_out(cacheline_t<W> * memptr, mask_t<W> & mask) {
    static_assert(false, "Not implemented");
  }
  
  /**
   * Brings consecutive words from the global memory starting at address 
   * \c memptr into this cacheline. \c memptr must be cacheline aligned.
   */
  void stream_in( W const * memptr) {
    static_assert(false, "Not implemented");
  }
  
  void stream_in(cacheline_t<W> const * memptr) {
    static_assert(false, "Not implemented");
  }
  
  /**
   * Sets all the words in this cacheline to zero.
   */
  void setzero() {
    static_assert(false, "Not implemented");
  }
  
  /**
   * Ensures that the content of this cacheline written using the stream_out() 
   * or masked_stream_out() operation is actually written to the global memory.
   */
  static void flush() {
    static_assert(false, "Not implemented");
  }
  
  /**
   * Returns a largest word pointer smaller than \c ptr aligned at a 
   * cacheline_t<W> boundary.
   */
  template<typename T>
  W* aligned_pointer(T *ptr) {
    static_assert(false, "Not implemented");
  }
  
  /**
   * Returns a largest \c const word pointer smaller than \c ptr aligned at a 
   * cacheline_t<W> boundary.
   */
  template<typename T>
  const W* aligned_pointer(const T *ptr) {
    static_assert(false, "Not implemented");
  }
  
  /**
   * Returns the number of words by which \c ptr is ahead of the largest 
   * cacheline_t<W> aligned pointer smaller than \c ptr.
   */
  template<typename PT>
  uintptr_t aligned_offset(PT ptr) {
    static_assert(false, "Not implemented");
  }
};

namespace cacheline_internal {

/**
 * Returns a largest word pointer smaller than \c ptr aligned at a 
 * cacheline_t<W> boundary.
 */
//template<typename PT, typename W>
//typename make_word_pointer<PT, W>::type aligned_pointer(PT ptr) {
//  return reinterpret_cast<typename make_word_pointer<PT, W>::type>(
//    ((uintptr_t)(void*)ptr & (uintptr_t)(-alignof(cacheline_t<W>)))
//  );
//}
template<typename PT, typename W>
typename convert_pointer<PT, cacheline_t<W>>::type
aligned_pointer(PT ptr) {
  return reinterpret_cast<
    typename convert_pointer<PT, cacheline_t<W>>::type
  >(
    qsm::aligned_pointer((void*)ptr, alignof(cacheline_t<W>))
  );
}

/**
 * Returns the number of words by which \c ptr is ahead of the largest 
 * cacheline_t<W> aligned pointer smaller than \c ptr.
 */
//template<typename PT, typename W>
//uintptr_t aligned_offset(PT ptr) {
//  return ((uintptr_t)(void*)ptr & (uintptr_t)(alignof(cacheline_t<W>)-1))/sizeof(W);
//}
template<typename PT, typename W>
uintptr_t aligned_offset(PT ptr) {
  return qsm::aligned_offset((void*)ptr, alignof(cacheline_t<W>))/sizeof(W);
}

// NOTE: Choosing uint8_t, uint16_t or uint32_t as word_t did not 
// show any performance difference on IB. On KNC, uint32_t is the 
// only word size since the instructions do not support manipulation 
// of char or short data.
/** Default type for a cacheline word. */
typedef uint8_t default_cacheline_word_type;
//typedef uint32_t word_t;
//typedef CL_WORD_TYPE word_t;

} // end: namespace cacheline_internal

// Include platform specific implementations
#if defined(__AVX__)
#include "cacheline_avx.hpp"
#else
#error "Platform specific implementation of cacheline.hpp not found"
#endif

/** Default type for a mask. */
typedef mask_t<cacheline_internal::default_cacheline_word_type> mask;

/** Default type for a cacheline. */
typedef cacheline_t<cacheline_internal::default_cacheline_word_type> cacheline;

/** Default type for a pointer to an instance of cacheline. */
//typedef cacheline_t<word_t>* __restrict__ __attribute__((align_value(64))) cacheline_ptr;

} // end: namespace qsm
