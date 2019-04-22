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

#include <cacheline.hpp>

#include <tuple>
#include <type_traits>

namespace qsm {

// GI - index type for the global array
// LI - index type for the look-back
// T  - type of the global array
template<
  typename GI,
  typename LI,
  typename... T
>
class cgrqueue {
private:
  std::tuple<T const * ...> m_gptr;
  GI m_gidx;
  
public:
  cgrqueue() : m_gptr{}, m_gidx(0) {
  }
  
  void configure(T const * ... gptr) {
    m_gptr = { gptr... };
    m_gidx = 0;
  }
  
  void configure(GI const gidx, T const * ... gptr) {
    m_gptr = { gptr... };
    m_gidx = gidx;
  }
  
  void configure(std::tuple<T const * ...> const & gptr) {
    m_gptr = gptr;
    m_gidx = 0;
  }
  
  void configure(GI const gidx, std::tuple<T const * ...> const & gptr) {
    m_gptr = gptr;
    m_gidx = gidx;
  }
  
  GI l2gindex(LI const lidx) {
    return m_gidx - lidx;
  }
  
  template<size_t I>
  typename std::tuple_element<I, std::tuple<T const & ...>>::type
  get(GI const idx) const {
    return std::get<I>(m_gptr)[idx];
  }
  
  void advance(GI const n) {
    m_gidx += n;
  }
  
  template<int I>
  char const * get_current_char_pointer() const {
    return reinterpret_cast<char const *>(&(std::get<I>(m_gptr)[m_gidx]));
  }
};


// GI - index type for the global array
// LI - index type for the look-back
// T  - type of the global array
template<
  typename GI,
  typename LI,
  typename T
>
class cgrqueue<GI, LI, T> {
private:
  T const * m_gptr;
  GI m_gidx;
  
public:
  cgrqueue() : m_gptr(nullptr), m_gidx(0) {
  }
  
  void configure(T const * gptr) {
    m_gptr = gptr;
    m_gidx = 0;
  }
  
  void configure(GI const gidx, T const * gptr) {
    m_gptr = gptr;
    m_gidx = gidx;
  }
  
  GI l2gindex(LI const lidx) {
    return m_gidx - lidx;
  }
  
  T const & get(GI const idx) const {
    return m_gptr[idx];
  }
  
  void advance(GI const n) {
    m_gidx += n;
  }
  
  char const * get_current_char_pointer() const {
    return reinterpret_cast<char const *>(&m_gptr[m_gidx]);
  }
};

// GI - index type for the global array
// T  - type of the global array
template<typename GI, typename... T>
class seq_cgrqueue {
private:
  std::tuple<T const * ...> m_gptr;
  GI m_gidx;
  GI m_mark;
  
public:
  seq_cgrqueue() : m_gptr{}, m_gidx(0), m_mark(0) {
  }
  
  void configure(T const * ... gptr) {
    m_gptr = { gptr... };
    m_gidx = 0;
    m_mark = 0;
  }
  
  void configure(GI const gidx, T const * ... gptr) {
    m_gptr = { gptr... };
    m_gidx = gidx;
    m_mark = gidx;
  }
  
  void configure(std::tuple<T const * ...> const & gptr) {
    m_gptr = gptr;
    m_gidx = 0;
    m_mark = 0;
  }
  
  void configure(GI const gidx, std::tuple<T const * ...> const & gptr) {
    m_gptr = gptr;
    m_gidx = gidx;
    m_mark = gidx;
  }
  
  template<size_t I>
  typename std::tuple_element<I, std::tuple<T const & ...>>::type
  get() const {
    return std::get<I>(m_gptr)[m_gidx];
  }
  
  void get(std::tuple<T * ...> const & op, GI idx) const {
    #pragma forceinline
    qsm::assign(op, idx, m_gptr, m_gidx);
  }
  
  //std::tuple<T const & ...> pop() {
  //  GI idx = m_gidx++;
  //  return qsm::apply_with_tuple_result(
  //    [idx](auto gptr) {
  //      return gptr[idx];
  //    }, m_gptr
  //  );
  //}
  
  void advance() {
    m_gidx++;
  }
  
  void mark() {
    m_mark = m_gidx;
  }
  
  void rewind() {
    m_gidx = m_mark;
  }
};

// GI - index type for the global array
// T  - type of the global array
template<typename GI, typename T>
class seq_cgrqueue<GI, T> {
private:
  T const * m_gptr;
  GI m_gidx;
  GI m_mark;
  
public:
  seq_cgrqueue() : m_gptr(nullptr), m_gidx(0), m_mark(0) {
  }
  
  void configure(T const * gptr) {
    m_gptr = gptr;
    m_gidx = 0;
    m_mark = 0;
  }
  
  void configure(GI const gidx, T const * gptr) {
    m_gptr = gptr;
    m_gidx = gidx;
    m_mark = gidx;
  }
  
  T const & get() const {
    return m_gptr[m_gidx];
  }
  
  T const & pop() {
    return m_gptr[m_gidx++];
  }
  
  void advance() {
    m_gidx++;
  }
  
  void mark() {
    m_mark = m_gidx;
  }
  
  void rewind() {
    m_gidx = m_mark;
  }
};

//////////////////////
// Reverse read queues
//////////////////////

// GI - index type for the global array
// LI - index type for the look-back
// T  - type of the global array
template<
  typename GI,
  typename LI,
  typename... T
>
class rcgrqueue {
private:
  std::tuple<T const * ...> m_gptr;
  GI m_offset;
  
public:
  rcgrqueue() : m_gptr{}, m_offset(0) {
  }
  
  void configure(T const * ... gptr) {
    m_gptr = { gptr... };
    m_offset = 0;
  }
  
  void configure(GI const offset, T const * ... gptr) {
    m_gptr = { gptr... };
    m_offset = offset;
  }
  
  void configure(std::tuple<T const * ...> const & gptr) {
    m_gptr = gptr;
    m_offset = 0;
  }
  
  void configure(GI const offset, std::tuple<T const * ...> const & gptr) {
    m_gptr = gptr;
    m_offset = offset;
  }
  
  GI l2gindex(LI const lidx) {
    return m_offset + lidx;
  }
  
  template<size_t I>
  typename std::tuple_element<I, std::tuple<T const & ...>>::type
  get(GI const idx) const {
    return std::get<I>(m_gptr)[idx];
  }
  
  void advance(GI const n) {
    m_offset -= n;
  }
  
  template<int I>
  char const * get_current_char_pointer() const {
    return reinterpret_cast<char const *>(&(std::get<I>(m_gptr)[m_offset]));
  }
};


// GI - index type for the global array
// LI - index type for the look-back
// T  - type of the global array
template<
  typename GI,
  typename LI,
  typename T
>
class rcgrqueue<GI, LI, T> {
private:
  T const * m_gptr;
  GI m_offset;
  
public:
  rcgrqueue() : m_gptr(nullptr), m_offset(0) {
  }
  
  void configure(T const * gptr) {
    m_gptr = gptr;
    m_offset = 0;
  }
  
  void configure(GI const offset, T const * gptr) {
    m_gptr = gptr;
    m_offset = offset;
  }
  
  GI l2gindex(LI const lidx) {
    return m_offset + lidx;
  }
  
  T const & get(GI const idx) const {
    return m_gptr[idx];
  }
  
  void advance(GI const n) {
    m_offset -= n;
  }
  
  char const * get_current_char_pointer() const {
    return reinterpret_cast<char const *>(&m_gptr[m_offset]);
  }
};


// GI - index type for the global array
// T  - type of the global array
template<typename GI, typename... T>
class seq_rcgrqueue {
private:
  std::tuple<T const * ...> m_gptr;
  GI m_offset;
  GI m_offset_mark;
  
public:
  seq_rcgrqueue() : m_gptr{}, m_offset(0), m_offset_mark(0) {
  }
  
  void configure(T const * ... gptr) {
    m_gptr = { gptr... };
    m_offset = 0;
    m_offset_mark = 0;
  }
  
  void configure(GI const offset, T const * ... gptr) {
    m_gptr = { gptr... };
    m_offset = offset;
    m_offset_mark = offset;
  }
  
  void configure(std::tuple<T const * ...> const & gptr) {
    m_gptr = gptr;
    m_offset = 0;
    m_offset_mark = 0;
  }
  
  void configure(GI const gidx, std::tuple<T const * ...> const & gptr) {
    m_gptr = gptr;
    m_offset = offset;
    m_offset_mark = offset;
  }
  
  template<size_t I>
  typename std::tuple_element<I, std::tuple<T const & ...>>::type
  get() const {
    return std::get<I>(m_gptr)[m_offset];
  }
  
  std::tuple<T const & ...> pop() {
    GI off = --m_offset;
    return qsm::apply_with_tuple_result(
      [off](auto gptr) {
        return gptr[off];
      }, m_gptr
    );
  }
  
  void advance() {
    --m_offset;
  }
  
  void mark() {
    m_offset_mark = m_offset;
  }
  
  void rewind() {
    m_offset = m_offset_mark;
  }
};


// GI - index type for the global array
// T  - type of the global array
template<typename GI, typename T>
class seq_rcgrqueue<GI, T> {
private:
  T const * m_gptr;
  GI m_offset;
  GI m_offset_mark;
  
public:
  seq_rcgrqueue() : m_gptr(nullptr), m_offset(0), m_offset_mark(0) {
  }
  
  void configure(T const * gptr) {
    m_gptr = gptr;
    m_offset = 0;
    m_offset_mark = 0;
  }
  
  void configure(GI const offset, T const * gptr) {
    m_gptr = gptr;
    m_offset = offset;
    m_offset_mark = offset;
  }
  
  T const & get() const {
    return m_gptr[m_offset];
  }
  
  T const & pop() {
    return m_gptr[--m_offset];
  }
  
  void advance() {
    --m_offset;
  }
  
  void mark() {
    m_offset_mark = m_offset;
  }
  
  void rewind() {
    m_offset = m_offset_mark;
  }
};

} // end: namespace qsm
