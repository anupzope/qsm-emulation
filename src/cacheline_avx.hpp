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

#include <cstdint>
#include <immintrin.h>

template<>
class mask_t<uint8_t> {
  friend class cacheline_t<uint8_t>;
  
public:
  mask_t() {
    reset();
  }
  
  mask_t(int step) {
    stepup(step);
  }
  
  void set() {
    __m128i zero = _mm_setzero_si128();
    __m128i ones = _mm_cmpeq_epi32(zero, zero);
    r128[0] = ones; r128[1] = ones; r128[2] = ones; r128[3] = ones;
  }
  
  void reset() {
    r256[0] = _mm256_setzero_si256();
    r256[1] = _mm256_setzero_si256();
  }
  
  void stepup(int step) {
    uint64_t mask = ((uint64_t)-1 << step);
    set(mask);
  }
  
  void stepdown(int step) {
    uint64_t mask = ((uint64_t)0x1u << step) - 1;
    set(mask);
  }
  
  void flip() {
    __m128i zero = _mm_setzero_si128();
    __m128i ones = _mm_cmpeq_epi8(zero, zero);
    r128[0] = _mm_andnot_si128(r128[0], ones);
    r128[1] = _mm_andnot_si128(r128[1], ones);
    r128[2] = _mm_andnot_si128(r128[2], ones);
    r128[3] = _mm_andnot_si128(r128[3], ones);
  }
  
  bool is_zero() const {
    return _mm256_testz_si256(r256[0], r256[0]) && 
      _mm256_testz_si256(r256[1], r256[1]);
  }
  
  mask_t<uint8_t> & operator^=(mask_t<uint8_t> const & rhs) {
    r128[0] = _mm_xor_si128(r128[0], rhs.r128[0]);
    r128[1] = _mm_xor_si128(r128[1], rhs.r128[1]);
    r128[2] = _mm_xor_si128(r128[2], rhs.r128[2]);
    r128[3] = _mm_xor_si128(r128[3], rhs.r128[3]);
    return *this;
  }
  
  mask_t<uint8_t> & operator&=(mask_t<uint8_t> const & rhs) {
    r128[0] = _mm_and_si128(r128[0], rhs.r128[0]);
    r128[1] = _mm_and_si128(r128[1], rhs.r128[1]);
    r128[2] = _mm_and_si128(r128[2], rhs.r128[2]);
    r128[3] = _mm_and_si128(r128[3], rhs.r128[3]);
    return *this;
  }
  
private:
  void set(uint64_t value) {
    __m128i m1 = _mm_set_epi8(0x80, 0x80, 0x40, 0x40, 0x20, 0x20, 0x10, 0x10, 
                              0x08, 0x08, 0x04, 0x04, 0x02, 0x02, 0x01, 0x01);
    __m128i m2 = _mm_set_epi8(15, 13, 11, 9, 7, 5, 3, 1, 
                              14, 12, 10, 8, 6, 4, 2, 0);
    for(int i = 0; i < 4; ++i) {
      __m128i a = _mm_set1_epi16(value >> 16*i);
      r128[i] = _mm_shuffle_epi8(_mm_cmpeq_epi8(_mm_and_si128(a, m1), m1), m2);
    }
  }
  
private:
  union {
    __m128i r128[4];
    __m256i r256[2];
  };
};

template<>
class alignas(CLSize) cacheline_t<uint8_t> {
public:
  typedef uint8_t word_type;
  static constexpr unsigned int ByteSize = CLSize;
  static constexpr unsigned int WordSize = ByteSize/sizeof(word_type);
  
  template<typename T = word_type>
  static T * round_down(void * ptr) {
    return reinterpret_cast<T *>(
      (uintptr_t)ptr & ~(uintptr_t)(ByteSize-1)
    );
  }
  
  template<typename T = word_type>
  static T * round_up(void * ptr) {
    return reinterpret_cast<T *>(
      ((uintptr_t)ptr + ByteSize) & ~(uintptr_t)(ByteSize-1)
    );
  }
  
public:
  // memptr must be aligned
  inline void stream_out(word_type * memptr) {
    __assume_aligned(memptr, CLSize);
    _mm256_stream_si256(reinterpret_cast<__m256i *>(memptr), r256[0]);
    _mm256_stream_si256(reinterpret_cast<__m256i *>(memptr+32), r256[1]);
  }
  
  inline void stream_out(cacheline_t<word_type> * memptr) {
    stream_out(reinterpret_cast<word_type *>(memptr));
  }
  
  // memptr must be aligned
  void masked_stream_out(word_type * memptr, mask_t<word_type> & mask) {
    __assume_aligned(memptr, CLSize);
    _mm_maskmoveu_si128(_mm256_castsi256_si128(r256[0]), mask.r128[0], reinterpret_cast<char*>(memptr));
    _mm_maskmoveu_si128(_mm256_extractf128_si256(r256[0], 1), mask.r128[1], reinterpret_cast<char*>(memptr)+16);
    _mm_maskmoveu_si128(_mm256_castsi256_si128(r256[1]), mask.r128[2], reinterpret_cast<char*>(memptr)+32);
    _mm_maskmoveu_si128(_mm256_extractf128_si256(r256[1], 1), mask.r128[3], reinterpret_cast<char*>(memptr)+48);
  }
  
  void masked_stream_out(cacheline_t<word_type> * memptr, mask_t<word_type> & mask) {
    masked_stream_out(reinterpret_cast<word_type *>(memptr), mask);
  }
  
  // memptr must be aligned
  void stream_in(word_type const * memptr) {
    __assume_aligned(memptr, CLSize);
    r256[0] = _mm256_load_si256(reinterpret_cast<__m256i const *>(memptr));
    r256[1] = _mm256_load_si256(reinterpret_cast<__m256i const *>(memptr+32));
  }
  
  void stream_in(cacheline_t<word_type> const * memptr) {
    stream_in(reinterpret_cast<word_type const *>(memptr));
  }
  
  void setzero() {
    r256[0] = _mm256_setzero_si256();
    r256[1] = _mm256_setzero_si256();
  }
  
  static void flush() {
    _mm_mfence();
  }
  
  template<typename PT>
  static typename convert_pointer<PT, cacheline_t<word_type>>::type
  aligned_pointer(PT ptr) {
    return cacheline_internal::aligned_pointer<PT, word_type>(ptr);
  }
  
  template<typename PT>
  static uintptr_t aligned_offset(PT ptr) {
    return cacheline_internal::aligned_offset<PT, word_type>(ptr);
  }
  
private:
  union {
    __m128i r128[4];
    __m256i r256[2];
  };
};

template<>
class mask_t<uint16_t> {
  friend class cacheline_t<uint16_t>;
  
public:
  mask_t() {
    reset();
  }
  
  mask_t(int step) {
    stepup(step);
  }
  
  void set() {
    __m128i zero = _mm_setzero_si128();
    __m128i ones = _mm_cmpeq_epi32(zero, zero);
    r128[0] = ones; r128[1] = ones; r128[2] = ones; r128[3] = ones;
  }
  
  void reset() {
    r256[0] = _mm256_setzero_si256();
    r256[1] = _mm256_setzero_si256();
  }
  
  void stepup(int step) {
    uint32_t mask = ((uint32_t)-1 << step);
    set(mask);
  }
  
  //movd mm0, eax
  //punpcklbw mm0, mm0
  //pshufw mm0, mm0, 0x00
  //pand mm0, [mask8040201008040201h]
  //pcmpeb mm0, [mask8040201008040201h]
  void stepdown(int step) {
    uint32_t mask = ((uint32_t)0x1u << step) - 1;
    set(mask);
  }
  
  void flip() {
    __m128i zero = _mm_setzero_si128();
    __m128i ones = _mm_cmpeq_epi32(zero, zero);
    r128[0] = _mm_andnot_si128(r128[0], ones);
    r128[1] = _mm_andnot_si128(r128[1], ones);
    r128[2] = _mm_andnot_si128(r128[2], ones);
    r128[3] = _mm_andnot_si128(r128[3], ones);
  }
  
  bool is_zero() const {
    return _mm256_testz_si256(r256[0], r256[0]) && 
      _mm256_testz_si256(r256[1], r256[1]);
  }
  
  mask_t<uint16_t> & operator^=(mask_t<uint16_t> const & rhs) {
    r128[0] = _mm_xor_si128(r128[0], rhs.r128[0]);
    r128[1] = _mm_xor_si128(r128[1], rhs.r128[1]);
    r128[2] = _mm_xor_si128(r128[2], rhs.r128[2]);
    r128[3] = _mm_xor_si128(r128[3], rhs.r128[3]);
    return *this;
  }
  
private:
  void set(uint32_t value) {
    __m128i m1 = _mm_set_epi16(0x0080, 0x0040, 0x0020, 0x0010, 0x0008, 0x0004, 0x0002, 0x0001);
    __m128i m2 = _mm_set_epi16(0x8000, 0x4000, 0x2000, 0x1000, 0x0800, 0x0400, 0x0200, 0x0100);
    r128[0] = _mm_cmpeq_epi16(_mm_and_si128(_mm_set1_epi16(value), m1), m1);
    r128[1] = _mm_cmpeq_epi16(_mm_and_si128(_mm_set1_epi16(value >> 8), m1), m1);
    r128[2] = _mm_cmpeq_epi16(_mm_and_si128(_mm_set1_epi16(value >> 16), m2), m2);
    r128[3] = _mm_cmpeq_epi16(_mm_and_si128(_mm_set1_epi16(value >> 24), m2), m2);
  }
  
private:
  union {
    __m128i r128[4];
    __m256i r256[2];
  };
};

template<>
class alignas(64) cacheline_t<uint16_t> {
public:
  typedef uint16_t word_type;
  static constexpr unsigned int ByteSize = 64;
  static constexpr unsigned int WordSize = 64/sizeof(word_type);
  
  template<typename T = word_type>
  static T * round_down(void * ptr) {
    return reinterpret_cast<T *>(
      (uintptr_t)ptr & ~(uintptr_t)(ByteSize-1)
    );
  }
  
  template<typename T = word_type>
  static T * round_up(void * ptr) {
    return reinterpret_cast<T *>(
      ((uintptr_t)ptr + ByteSize) & ~(uintptr_t)(ByteSize-1)
    );
  }
  
public:
  // memptr must be aligned
  void stream_out(word_type * memptr) {
    __assume_aligned(memptr, 64);
    _mm256_stream_si256(reinterpret_cast<__m256i*>(memptr), r256[0]);
    _mm256_stream_si256(reinterpret_cast<__m256i*>(memptr+16), r256[1]);
  }
  
  void stream_out(cacheline_t<word_type> * memptr) {
    stream_out(reinterpret_cast<word_type *>(memptr));
  }
  
  // memptr must be aligned
  void masked_stream_out(word_type * memptr, mask_t<word_type> & mask) {
    __assume_aligned(memptr, 64);
    _mm_maskmoveu_si128(_mm256_castsi256_si128(r256[0]), mask.r128[0], reinterpret_cast<char*>(memptr));
    _mm_maskmoveu_si128(_mm256_extractf128_si256(r256[0], 1), mask.r128[1], reinterpret_cast<char*>(memptr+8));
    _mm_maskmoveu_si128(_mm256_castsi256_si128(r256[1]), mask.r128[2], reinterpret_cast<char*>(memptr+16));
    _mm_maskmoveu_si128(_mm256_extractf128_si256(r256[1], 1), mask.r128[3], reinterpret_cast<char*>(memptr+24));
  }
  
  void masked_stream_out(cacheline_t<word_type> * memptr, mask_t<word_type> & mask) {
    masked_stream_out(reinterpret_cast<word_type *>(memptr), mask);
  }
  
  // memptr must be aligned
  void stream_in(word_type const * memptr) {
    __assume_aligned(memptr, 64);
    r256[0] = _mm256_load_si256(reinterpret_cast<__m256i const *>(memptr));
    r256[1] = _mm256_load_si256(reinterpret_cast<__m256i const *>(memptr+16));
  }
  
  void stream_in(cacheline_t<word_type> const * memptr) {
    stream_in(reinterpret_cast<word_type const *>(memptr));
  }
  
  void setzero() {
    r256[0] = _mm256_setzero_si256();
    r256[1] = _mm256_setzero_si256();
  }
  
  static void flush() {
    _mm_mfence();
  }
  
  template<typename PT>
  static typename convert_pointer<PT, cacheline_t<word_type>>::type
  aligned_pointer(PT ptr) {
    return cacheline_internal::aligned_pointer<PT, word_type>(ptr);
  }
  
  template<typename PT>
  static uintptr_t aligned_offset(PT ptr) {
    return cacheline_internal::aligned_offset<PT, word_type>(ptr);
  }
  
private:
  union {
    __m128i r128[4];
    __m256i r256[2];
  };
};

template<>
class mask_t<uint32_t> {
  friend class cacheline_t<uint32_t>;
  
public:
  mask_t() {
    reset();
  }
  
  mask_t(int step) {
    stepup(step);
  }
  
  void set() {
    __m128i zero = _mm_setzero_si128();
    __m128i ones = _mm_cmpeq_epi32(zero, zero);
    r128[0] = ones; r128[1] = ones; r128[2] = ones; r128[3] = ones;
  }
  
  void reset() {
    r256[0] = _mm256_setzero_si256();
    r256[1] = _mm256_setzero_si256();
  }
  
  void stepup(int step) {
    uint16_t mask = ((uint16_t)-1 << step);
    set(mask);
  }
  
  void stepdown(int step) {
    uint16_t mask = ((uint16_t)0x1u << step) - 1;
    set(mask);
  }
  
  void flip() {
    __m128i ones = _mm_set1_epi32((uint32_t)-1);
    r128[0] = _mm_andnot_si128(r128[0], ones);
    r128[1] = _mm_andnot_si128(r128[1], ones);
    r128[2] = _mm_andnot_si128(r128[2], ones);
    r128[3] = _mm_andnot_si128(r128[3], ones);
  }
  
  bool is_zero() const {
    return _mm256_testz_si256(r256[0], r256[0]) && 
      _mm256_testz_si256(r256[1], r256[1]);
  }
  
  mask_t<uint32_t> & operator^=(mask_t<uint32_t> const & rhs) {
    r128[0] = _mm_xor_si128(r128[0], rhs.r128[0]);
    r128[1] = _mm_xor_si128(r128[1], rhs.r128[1]);
    r128[2] = _mm_xor_si128(r128[2], rhs.r128[2]);
    r128[3] = _mm_xor_si128(r128[3], rhs.r128[3]);
    return *this;
  }
  
private:
  void set(uint16_t value) {
    __m128i op = _mm_set1_epi32(value);
    __m128i m1 = _mm_set_epi32(0x8, 0x4, 0x2, 0x1);
    __m128i m2 = _mm_set_epi32(0x80, 0x40, 0x20, 0x10);
    __m128i m3 = _mm_set_epi32(0x800, 0x400, 0x200, 0x100);
    __m128i m4 = _mm_set_epi32(0x8000, 0x4000, 0x2000, 0x1000);
    r128[0] = _mm_cmpeq_epi16(_mm_and_si128(op, m1), m1);
    r128[1] = _mm_cmpeq_epi16(_mm_and_si128(op, m2), m2);
    r128[2] = _mm_cmpeq_epi16(_mm_and_si128(op, m3), m3);
    r128[3] = _mm_cmpeq_epi16(_mm_and_si128(op, m4), m4);
  }
  
private:
  union {
    __m128i r128[4];
    __m256i r256[2];
  };
};

template<>
class alignas(64) cacheline_t<uint32_t> {
public:
  typedef uint32_t word_type;
  static constexpr unsigned int ByteSize = 64;
  static constexpr unsigned int WordSize = 64/sizeof(word_type);
  
  template<typename T = word_type>
  static T * round_down(void * ptr) {
    return reinterpret_cast<T *>(
      (uintptr_t)ptr & ~(uintptr_t)(ByteSize-1)
    );
  }
  
  template<typename T = word_type>
  static T * round_up(void * ptr) {
    return reinterpret_cast<T *>(
      ((uintptr_t)ptr + ByteSize) & ~(uintptr_t)(ByteSize-1)
    );
  }
  
public:
  // memptr must be aligned
  void stream_out(word_type * memptr) {
    __assume_aligned(memptr, 64);
    _mm256_stream_si256(reinterpret_cast<__m256i*>(memptr), r256[0]);
    _mm256_stream_si256(reinterpret_cast<__m256i*>(memptr+8), r256[1]);
  }
  
  void stream_out(cacheline_t<word_type> * memptr) {
    stream_out(reinterpret_cast<word_type *>(memptr));
  }
  
  // memptr must be aligned
  void masked_stream_out(word_type * memptr, mask_t<word_type> &mask) {
    __assume_aligned(memptr, 64);
    _mm_maskmoveu_si128(_mm256_castsi256_si128(r256[0]), mask.r128[0], reinterpret_cast<char*>(memptr));
    _mm_maskmoveu_si128(_mm256_extractf128_si256(r256[0], 1), mask.r128[1], reinterpret_cast<char*>(memptr+4));
    _mm_maskmoveu_si128(_mm256_castsi256_si128(r256[1]), mask.r128[2], reinterpret_cast<char*>(memptr+8));
    _mm_maskmoveu_si128(_mm256_extractf128_si256(r256[1], 1), mask.r128[3], reinterpret_cast<char*>(memptr+12));
  }
  
  void masked_stream_out(cacheline_t<word_type> * memptr, mask_t<word_type> &mask) {
    masked_stream_out(reinterpret_cast<word_type *>(memptr), mask);
  }
  
  // memptr must be aligned
  void stream_in(word_type const * memptr) {
    __assume_aligned(memptr, 64);
    r256[0] = _mm256_load_si256(reinterpret_cast<__m256i const *>(memptr));
    r256[1] = _mm256_load_si256(reinterpret_cast<__m256i const *>(memptr+8));
  }
  
  void stream_in(cacheline_t<word_type> const * memptr) {
    stream_in(reinterpret_cast<word_type const *>(memptr));
  }
  
  void setzero() {
    r256[0] = _mm256_setzero_si256();
    r256[1] = _mm256_setzero_si256();
  }
  
  static void flush() {
    _mm_mfence();
  }
  
  template<typename PT>
  static typename convert_pointer<PT, cacheline_t<word_type>>::type
  aligned_pointer(PT ptr) {
    return cacheline_internal::aligned_pointer<PT, word_type>(ptr);
  }
  
  template<typename PT>
  static uintptr_t aligned_offset(PT ptr) {
    return cacheline_internal::aligned_offset<PT, word_type>(ptr);
  }
  
private:
  union {
    __m128i r128[4];
    __m256i r256[2];
  };
};
