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

#include "debruijn.hpp"
#include <immintrin.h>

namespace qsm {

template<>
inline uint8_t bit1atn(int index) {
  return (uint8_t)0x1u << index;
}

template<>
inline uint16_t bit1atn(int index) {
  return (uint16_t)0x1u << index;
}

template<>
inline uint32_t bit1atn(int index) {
  return (uint32_t)0x1u << index;
}

template<>
inline uint64_t bit1atn(int index) {
  return (uint64_t)0x1u << index;
}

template<>
inline __m128i bit1atn<__m128i>(int index) {
  __m128i uf = _mm_set_epi32(128, 96, 64, 32);
  __m128i lf = _mm_set_epi32(95, 63, 31, -1);
  
  __m128i idx = _mm_set1_epi32(index);
  
  __m128i zero = _mm_setzero_si128();
  
  __m128i allone;
  allone = _mm_cmpeq_epi32(allone, allone);
  
  __m128i one = _mm_srli_epi32(allone, 31);
  
  __m128i shiftmask = _mm_srli_epi64(allone, 59);
  __m128i shiftvalue = _mm_sll_epi32(one, _mm_and_si128(idx, shiftmask));
  
  __m128i blendmask = _mm_and_si128(_mm_cmplt_epi32(idx, uf), _mm_cmpgt_epi32(idx, lf));
  return _mm_blendv_epi8(zero, shiftvalue, blendmask);
}

template<>
inline __m256i bit1atn<__m256i>(int index) {
  return _mm256_set_m128i(bit1atn<__m128i>(index-128), bit1atn<__m128i>(index));
}




template<>
inline uint8_t setbit(uint8_t v, int index) {
  return v | bit1atn<uint8_t>(index);
}

template<>
inline uint16_t setbit(uint16_t v, int index) {
  return v | bit1atn<uint16_t>(index);
}

template<>
inline uint32_t setbit(uint32_t v, int index) {
  return v | bit1atn<uint32_t>(index);
}

template<>
inline uint64_t setbit(uint64_t v, int index) {
  return v | bit1atn<uint64_t>(index);
}

template<>
inline __m128i setbit<__m128i>(__m128i v, int index) {
  return _mm_or_si128(v, bit1atn<__m128i>(index));
}

template<>
inline __m256i setbit<__m256i>(__m256i v, int index) {
  __m128i lo = setbit<__m128i>(_mm256_castsi256_si128(v), index);
  __m128i hi = setbit<__m128i>(_mm256_castsi256_si128(_mm256_permute2f128_si256(v, v, 0x11)), index-128);
  return _mm256_set_m128i(hi, lo);
}




template<>
inline uint8_t setzero() {
  return 0;
}

template<>
inline uint16_t setzero() {
  return 0;
}
template<>
inline uint32_t setzero() {
  return 0;
}
template<>
inline uint64_t setzero() {
  return 0;
}

template<>
inline __m128i setzero<__m128i>() {
  return _mm_setzero_si128();
}

template<>
inline __m256i setzero<__m256i>() {
  return _mm256_setzero_si256();
}




template<>
inline bool testzero(uint8_t v) {
  return !v;
}

template<>
inline bool testzero(uint16_t v) {
  return !v;
}

template<>
inline bool testzero(uint32_t v) {
  return !v;
}

template<>
inline bool testzero(uint64_t v) {
  return !v;
}

template<>
inline bool testzero<__m128i>(__m128i v) {
  return _mm_testz_si128(v, v);
}

template<>
inline bool testzero<__m256i>(__m256i v) {
  return _mm256_testz_si256(v, v);
}




template<>
inline bool testnotzero(uint8_t v) {
  return v;
}

template<>
inline bool testnotzero(uint16_t v) {
  return v;
}

template<>
inline bool testnotzero(uint32_t v) {
  return v;
}

template<>
inline bool testnotzero(uint64_t v) {
  return v;
}

template<>
inline bool testnotzero<__m128i>(__m128i v) {
  return !_mm_testz_si128(v, v);
}

template<>
inline bool testnotzero<__m256i>(__m256i v) {
  return !_mm256_testz_si256(v, v);
}




template<>
inline bool testnotequal(uint8_t a, uint8_t b) {
  return a != b;
}

template<>
inline bool testnotequal(uint16_t a, uint16_t b) {
  return a != b;
}

template<>
inline bool testnotequal(uint32_t a, uint32_t b) {
  return a != b;
}

template<>
inline bool testnotequal(uint64_t a, uint64_t b) {
  return a != b;
}

template<>
inline bool testnotequal<__m128i>(__m128i a, __m128i b) {
  return _mm_movemask_epi8(_mm_cmpeq_epi8(a, b)) != 0xFFFF;
}

template<>
inline bool testnotequal<__m256i>(__m256i a, __m256i b) {
  return _mm256_movemask_ps(_mm256_cmp_ps(_mm256_castsi256_ps(a), _mm256_castsi256_ps(b), _CMP_EQ_US)) != 0xFF;
}




template<>
inline bool testequal(uint8_t a, uint8_t b) {
  return a == b;
}

template<>
inline bool testequal(uint16_t a, uint16_t b) {
  return a == b;
}

template<>
inline bool testequal(uint32_t a, uint32_t b) {
  return a == b;
}

template<>
inline bool testequal(uint64_t a, uint64_t b) {
  return a == b;
}

template<>
inline bool testequal<__m128i>(__m128i a, __m128i b) {
  return _mm_movemask_epi8(_mm_cmpeq_epi8(a, b)) == 0xFFFF;
}

template<>
inline bool testequal<__m256i>(__m256i a, __m256i b) {
  return _mm256_movemask_ps(_mm256_cmp_ps(_mm256_castsi256_ps(a), _mm256_castsi256_ps(b), _CMP_EQ_US)) == 0xFF;
}




template<>
inline int bit1count<uint8_t>(uint8_t v) {
#ifdef __KNC__
  return _mm_countbits_32(v);
#else
  return _mm_popcnt_u32(v);
#endif
}

template<>
inline int bit1count<uint16_t>(uint16_t v) {
#ifdef __KNC__
  return _mm_countbits_32(v);
#else
  return _mm_popcnt_u32(v);
#endif
}

template<>
inline int bit1count<uint32_t>(uint32_t v) {
#ifdef __KNC__
  return _mm_countbits_32(v);
#else
  return _mm_popcnt_u32(v);
#endif
}

template<>
inline int bit1count<uint64_t>(uint64_t v) {
#ifdef __KNC__
  return _mm_countbits_64(v);
#else
  return _mm_popcnt_u64(v);
#endif
}

template<>
inline int bit1count<__m128i>(__m128i v) {
  const __m128i vhi = _mm_unpackhi_epi64(v, v);
  return _popcnt64(_mm_cvtsi128_si64(v)) + _popcnt64(_mm_cvtsi128_si64(vhi));
}

template<>
inline int bit1count<__m256i>(__m256i v) {
  __m256i vhi = _mm256_permute2f128_si256(v, v, 0x11);
  return bit1count(_mm256_castsi256_si128(v)) + bit1count(_mm256_castsi256_si128(vhi));
}




#define BSC_LS1BINDEX 0

template<>
inline int ls1bindex<uint8_t>(uint8_t v) {
#if BSC_LS1BINDEX == 0
  return _bit_scan_forward(v);
#elif BSC_LS1BINDEX == 1
  return debruijn8table_2n[(uint8_t)((uint8_t)(v&-v) * debruijn8_2n) >> 5];
#else
  static_assert(false, "Not implemented")
#endif
}

template<>
inline int ls1bindex<uint16_t>(uint16_t v) {
#if BSC_LS1BINDEX == 0
  return _bit_scan_forward(v);
#elif BSC_LS1BINDEX == 1
  return debruijn16table_2n[(uint16_t)((uint16_t)(v&-v) * debruijn16_2n) >> 12];
#else
  static_assert(false, "Not implemented")
#endif
}

template<>
inline int ls1bindex<uint32_t>(uint32_t v) {
#if BSC_LS1BINDEX == 0
  return _bit_scan_forward(v);
#elif BSC_LS1BINDEX == 1
  return debruijn32table_2n[(uint32_t)((uint32_t)(v&-v) * debruijn32_2n) >> 27];
#else
  static_assert(false, "Not implemented");
#endif
}

template<>
inline int ls1bindex<uint64_t>(uint64_t v) {
#if BSC_LS1BINDEX == 0
  // TODO: remove these conditionals
  unsigned __int32 t;
  if(_BitScanForward64(&t, v)) {
    return t;
  }
  return 0;
#elif BSC_LS1BINDEX == 1
  return debruijn64table_2n[(uint64_t)((uint64_t)(v&-v) * debruijn64_2n) >> 58];
#else
  static_assert(false, "Not implemented");
#endif
}

//#include <iostream>
//void print(__m128i * v) {
//  uint64_t *p = reinterpret_cast<uint64_t*>(v);
//  std::cerr << "v = " << p[0] << ", " << p[1] << std::endl;
//}

template<>
inline int ls1bindex<__m128i>(__m128i v) {
  static int table[4] = {0, 0, 1, 0};
  __m128i zero = _mm_setzero_si128();
  __m128i ones = _mm_cmpeq_epi64(zero, zero);
  __m128i m = _mm_andnot_si128(_mm_cmpeq_epi64(v, zero), ones);
  int msk = _mm_movemask_epi8(m);
  int i = ((msk & 0x00FF) == 0x00FF) | (((msk & 0xFF00) == 0xFF00) << 1);
  int j = table[i];
  uint64_t *p = reinterpret_cast<uint64_t*>(&v);
  return debruijn64table_2n[(uint64_t)((uint64_t)(p[j]&-p[j]) * debruijn64_2n) >> 58];
}

/*template<>
inline int ls1bindex<__m256i>(__m256i v) {
  static int table[16] = {
  };
}*/

#define BSC_LOG2 0

template<>
inline uint8_t log2<uint8_t>(uint8_t v) {
#if BSC_LOG2 == 0
  return _bit_scan_reverse(v);
#elif BSC_LOG2 == 1
  uint8_t r = 0, s;
  s = (v > 0xF) << 2; v >>= r;
  s = (v > 0x3) << 1; v >>= s; r |= s;
                               r |= (v >> 1);
  return r;
#elif BSC_LOG2 == 2
  static_assert(false, "Not implemented");
  //v |= v >> 1; // first round down to one less than a power of 2 
  //v |= v >> 2;
  //v |= v >> 4;
  //return debruijn8table_2nm1[(uint8_t)(v * debruijn8_2nm1) >> 5];
#else
  static_assert(false, "Not implemented");
#endif
}

template<>
inline uint16_t log2<uint16_t>(uint16_t v) {
#if BSC_LOG2 == 0
  return _bit_scan_reverse(v);
#elif BSC_LOG2 == 1
  uint16_t r = 0, s;
  s = (v > 0xFF) << 3; v >>= r;
  s = (v > 0xF ) << 2; v >>= r; r |= s;
  s = (v > 0x3 ) << 1; v >>= s; r |= s;
                                r |= (v >> 1);
  return r;
#elif BSC_LOG2 == 2
  static_assert(false, "Not implemented");
  //v |= v >> 1; // first round down to one less than a power of 2 
  //v |= v >> 2;
  //v |= v >> 4;
  //v |= v >> 8;
  //return debruijn16table_2nm1[(uint16_t)(v * debruijn16_2nm1) >> 12];
#else
  static_assert(false, "Not implemented");
#endif
}

template<>
inline uint32_t log2<uint32_t>(uint32_t v) {
#if BSC_LOG2 == 0
  return _bit_scan_reverse(v);
#elif BSC_LOG2 == 1
  uint32_t r = 0, s;
  r = (v > 0xFFFF) << 4; v >>= r;
  s = (v > 0xFF  ) << 3; v >>= s; r |= s;
  s = (v > 0xF   ) << 2; v >>= s; r |= s;
  s = (v > 0x3   ) << 1; v >>= s; r |= s;
                                  r |= (v >> 1);
  return r;
#elif BSC_LOG2 == 2
  v |= v >> 1; // first round down to one less than a power of 2 
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  return debruijn32table_2nm1[(uint32_t)(v * debruijn32_2nm1) >> 27];
#else
  static_assert(false, "Not implemented");
#endif
}

template<>
inline uint64_t log2<uint64_t>(uint64_t v) {
#if BSC_LOG2 == 0
  // TODO: remove these conditionals
  unsigned __int32 t;
  if(_BitScanReverse64(&t, v))
    return t;
  else
    return (uint64_t)0;
#elif BSC_LOG2 == 1
  uint64_t r = 0, s;
  r = (v > 0xFFFFFFFF) << 5; v >>= r;
  s = (v > 0xFFFF    ) << 4; v >>= s; r |= s;
  s = (v > 0xFF      ) << 3; v >>= s; r |= s;
  s = (v > 0xF       ) << 2; v >>= s; r |= s;
  s = (v > 0x3       ) << 1; v >>= s; r |= s;
                                      r |= (v >> 1);
  return r;
#elif BSC_LOG2 == 2
  v |= v >> 1; // first round down to one less than a power of 2 
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v |= v >> 32;
  return debruijn64table_2nm1[(uint64_t)(v * debruijn64_2nm1) >> 58];
#else
  static_assert(false, "Not implemented");
#endif
}

} // end: namespace qsm
