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

#include <array>

namespace qsm {

// GI - index type for the global array
// T  - type of the global array
template<typename GI, typename T, int BufferSize>
class cgwqueue {
  static_assert(sizeof(T)%sizeof(cacheline::word_type) == 0,
    "Invalid template arguments: size of T must be a multiple of size of word_t"
  );
  static_assert(alignof(T)%alignof(cacheline::word_type) == 0,
    "Invalid template arguments: alignment of T must be a multiple of alignment of word_t"
  );
  
private:
  static constexpr int BufferCachelineSize = integer_ceil(
    (size_t)BufferSize*sizeof(T), (size_t)cacheline::ByteSize
  ) + 1;
  
  static constexpr int BufferWordSize = BufferCachelineSize*cacheline::WordSize;
  
public:
  using write_pointer_type = T *;
  
private:
  using cacheline_pointer_type = typename convert_pointer<write_pointer_type, cacheline>::type;
  
  using word_pointer_type = typename convert_pointer<write_pointer_type, cacheline::word_type>::type;
  
private:
  cacheline_pointer_type m_gptr;
  GI m_gidx; // current index in global memory (in cachelines)
  int m_bidx; // current index in buffer (in words)
  cacheline m_buf[BufferCachelineSize];
  //cacheline_pointer_type m_buf;
  mask m_mask;
  
public:
  cgwqueue() : m_gptr(nullptr), m_gidx(0), m_bidx(0) {
  }
  
  void configure(write_pointer_type gptr/*, cacheline_pointer_type cptr*/) {
    m_gptr = cacheline::aligned_pointer(gptr);
    //m_buf = cptr;
    m_gidx = 0;
    m_bidx = cacheline::aligned_offset(gptr);
    m_mask.stepup(m_bidx);
  }
  
  int get_buffer_size() const {
    return BufferSize;
  }
  
  write_pointer_type get(int * bsz) {
    *bsz = BufferSize;
    return reinterpret_cast<write_pointer_type>(
      &(reinterpret_cast<word_pointer_type>(m_buf)[m_bidx])
    );
  }
  
  write_pointer_type get() {
    return reinterpret_cast<write_pointer_type>(
      &(reinterpret_cast<word_pointer_type>(m_buf)[m_bidx])
    );
  }
  
  void flush(int const count) {
    int const nwords = count*sizeof(T)/sizeof(cacheline::word_type);
    int const ebidx = m_bidx + nwords;
    int const ebclidx = ebidx/cacheline::WordSize;
    int const ebcloff = ebidx%cacheline::WordSize;
    int clidx = 0;
    
    if(!m_gidx && clidx < ebclidx) {
      m_buf[clidx++].masked_stream_out(&m_gptr[m_gidx++], m_mask);
      m_mask.set();
    }
    
    while(clidx < ebclidx) {
      m_buf[clidx++].stream_out(&m_gptr[m_gidx++]);
    }
    
    m_buf[0] = m_buf[ebclidx];
    m_bidx = ebcloff;
  }
  
  void deactivate() {
    if(m_bidx) {
      m_mask ^= mask(m_bidx);
      m_buf[0].masked_stream_out(&m_gptr[m_gidx++], m_mask);
      
      m_gptr = nullptr;
      m_gidx = 0;
      m_bidx = 0;
      m_mask.reset();
    }
  }
};

// GI - index type for the global array
// T  - type of the global array
template<typename GI, typename T, int BufferSize>
class rcgwqueue {
  static_assert(sizeof(T)%sizeof(cacheline::word_type) == 0,
    "Invalid template arguments: size of T must be a multiple of size of word_t"
  );
  static_assert(alignof(T)%alignof(cacheline::word_type) == 0,
    "Invalid template arguments: alignment of T must be a multiple of alignment of word_t"
  );
  
private:
  static constexpr int BufferCachelineSize = integer_ceil(
    (size_t)BufferSize*sizeof(T), (size_t)cacheline::ByteSize
  ) + 1;
  
  static constexpr int BufferWordSize = BufferCachelineSize*cacheline::WordSize;
  
public:
  using write_pointer_type = T *;
  
private:
  using cacheline_pointer_type = typename convert_pointer<write_pointer_type, cacheline>::type;
  
  using word_pointer_type = typename convert_pointer<write_pointer_type, cacheline::word_type>::type;
  
private:
  cacheline_pointer_type m_gptr;
  GI m_gidx; // current index in global memory (in cachelines)
  int m_bidx; // current index in buffer (in words)
  cacheline m_buf[BufferCachelineSize];
  mask m_mask;
  
public:
  rcgwqueue() : m_gptr(nullptr), m_gidx(0), m_bidx(0) {
  }
  
  void configure(write_pointer_type gptr) {
    m_gidx = 0;
    m_gptr = cacheline::aligned_pointer(gptr);
    m_bidx = (BufferCachelineSize-1)*cacheline::WordSize+cacheline::aligned_offset(gptr);
    m_mask.stepdown(cacheline::aligned_offset(gptr));
  }
  
  write_pointer_type get() {
    return reinterpret_cast<write_pointer_type>(
      &(reinterpret_cast<word_pointer_type>(m_buf)[m_bidx])
    );
  }
  
  void flush(int const count) {
    int const nwords = count*sizeof(T)/sizeof(cacheline::word_type);
    int const ebidx = m_bidx - nwords;
    int const ebclidx = ebidx/cacheline::WordSize;
    int const ebcloff = ebidx%cacheline::WordSize;
    int clidx = m_bidx/cacheline::WordSize;
    
    if(!m_gidx && clidx > ebclidx) {
      m_buf[clidx--].masked_stream_out(&m_gptr[-m_gidx++], m_mask);
      m_mask.set();
    }
    
    while(clidx > ebclidx) {
      m_buf[clidx--].stream_out(&m_gptr[-m_gidx++]);
    }
    
    m_buf[BufferCachelineSize-1] = m_buf[ebclidx];
    m_bidx = (BufferCachelineSize-1)*cacheline::WordSize+ebcloff;
  }
  
  void deactivate() {
    int clidx = m_bidx/cacheline::WordSize;
    int cloff = m_bidx%cacheline::WordSize;
    m_mask &= mask(cloff);
    m_buf[clidx].masked_stream_out(&m_gptr[-m_gidx], m_mask);
    
    m_gptr = nullptr;
    m_gidx = 0;
    m_bidx = 0;
    m_mask.reset();
  }
};







// GI - index type for the global array
// T  - type of the global array
template<int BufferSize, typename GI, typename... T>
class cgwqueue_pack {
  static_assert(
    qsm::all_true<(sizeof(T)%sizeof(cacheline::word_type) == 0)...>::value,
    "invalid template arguments: size of T must be a multiple of size of word_t"
  );
  
  static_assert(
    qsm::all_true<(alignof(T)%alignof(cacheline::word_type) == 0)...>::value,
    "invalid template arguments: alignment of T must be a multiple of alignment of word_t"
  );
  
private:
  
  GI m_gidx[sizeof...(T)];
  cacheline * m_gptr[sizeof...(T)];
  int m_bidx[sizeof...(T)];
  std::tuple<
    cacheline[
      integer_ceil((size_t)BufferSize*sizeof(T), (size_t)cacheline::ByteSize)+1
    ]...
  > m_buf;
  mask m_mask[sizeof...(T)];
  
  template<size_t... I>
  void configure_detail(
    std::tuple<T * ...> const gptr,
    std::index_sequence<I...>
  ) {
    auto t1 = { (m_gptr[I] = cacheline::aligned_pointer(std::get<I>(gptr))) ... };
    auto t2 = { (m_bidx[I] = cacheline::aligned_offset(std::get<I>(gptr))) ... };
    for(int i = 0; i < sizeof...(I); ++i) {
      m_gidx[i] = 0;
      m_mask[i].stepup(m_bidx[i]);
    };
  }
  
  template<size_t... I>
  std::tuple<T * ...> get_detail(std::index_sequence<I...>) {
    return std::make_tuple(
      reinterpret_cast<typename std::tuple_element<I, std::tuple<T...>>::type *>(
        reinterpret_cast<cacheline::word_type *>(std::get<I>(m_buf))+m_bidx[I]
      )...
    );
  }
  
  template<size_t I>
  int flush_detail_detail(int const count) {
    int const nwords = count*sizeof(
        typename std::tuple_element<I, std::tuple<T...>>::type
      )/sizeof(cacheline::word_type);
    int const ebidx = m_bidx[I] + nwords;
    int const ebclidx = ebidx/cacheline::WordSize;
    int const ebcloff = ebidx%cacheline::WordSize;
    int clidx = 0;
    
    if(!m_gidx[I] && clidx < ebclidx) {
      std::get<I>(m_buf)[clidx++].masked_stream_out(
       &m_gptr[I][m_gidx[I]++], m_mask[I]
      );
      m_mask[I].set();
    }
    
    while(clidx < ebclidx) {
      std::get<I>(m_buf)[clidx++].stream_out(&m_gptr[I][m_gidx[I]++]);
    }
    
    std::get<I>(m_buf)[0] = std::get<I>(m_buf)[ebclidx];
    m_bidx[I] = ebcloff;
    
    return 0;
  }
  
  template<size_t... I>
  void flush_detail(int const count, std::index_sequence<I...>) {
    auto t = {
      flush_detail_detail<I>(count) ...
    };
  }
  
  template<size_t I>
  int deactivate_detail_detail() {
    if(m_bidx[I]) {
      m_mask[I] ^= mask(m_bidx[I]);
      std::get<I>(m_buf)[0].masked_stream_out(
        &m_gptr[I][m_gidx[I]], m_mask[I]
      );
      m_gptr[I] = nullptr;
      m_gidx[I] = 0;
      m_bidx[I] = 0;
      m_mask[I].reset();
    }
  }
  
  template<size_t... I>
  void deactivate_detail(std::index_sequence<I...>) {
    auto t = {
      deactivate_detail_detail<I>() ...
    };
  }
  
public:
  void configure(std::tuple<T * ...> const & gptr) {
    configure_detail(gptr, std::make_index_sequence<sizeof...(T)>{});
  }
  
  std::tuple<T *...> get() {
    return get_detail(std::make_index_sequence<sizeof...(T)>{});
  }
  
  void flush(int const count) {
    flush_detail(count, std::make_index_sequence<sizeof...(T)>{});
  }
  
  void deactivate() {
    deactivate_detail(std::make_index_sequence<sizeof...(T)>{});
  }
};

} // end: namespace qsm
