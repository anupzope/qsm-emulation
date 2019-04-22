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

#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <type_traits>
#include <tuple>

template<typename GI>
struct ispace {
  GI start, end, size;
  
  ispace() : start(0), end(0), size(0) {
  }
  
  ispace(GI s, GI e) : start(s), end(e), size(e-s) {}
  
  void set(GI s, GI e) {
    start = s;
    end = e;
    size = e-s;
  }
};

namespace qsm {

template<typename T>
constexpr T integer_ceil(T const numerator, T const denominator) {
  return numerator/denominator + (numerator%denominator ? 1 : 0);
}



template<typename T>
constexpr T integer_log2(T v) {
  unsigned r = 0;
  while(v >>= 1) { r++; }
  return r;
}



template<typename T>
constexpr bool is_integer_power_of_two(T v) {
  return (v & (v - 1)) == 0;
}



template<typename T>
void indexspace_partition(T const start, T const end, T const nparts,
T const partidx, T * part_start, T * part_end) {
  T const length = end-start;
  T div = length/nparts;
  T rem = length%nparts;
  *part_start = start + partidx*div + (partidx < rem ? partidx : rem);
  *part_end = start + (partidx+1)*div + ((partidx+1) < rem ? (partidx+1) : rem);
}

void * aligned_pointer(void * ptr, size_t alignment) {
  return (void*)((uintptr_t)ptr & (uintptr_t)-alignment);
}

uintptr_t aligned_offset(void * ptr, size_t alignment) {
  return ((uintptr_t)ptr & (uintptr_t)(alignment-1));
}



template<typename SPT, typename DT>
struct convert_pointer {
};

template<typename ST, typename DT>
struct convert_pointer<ST const *, DT> {
  static_assert(!std::is_pointer<ST>::value, "ST cannot be a pointer type");
  static_assert(!std::is_pointer<DT>::value, "DT cannot be a pointer type");
  using type = DT const *;
};

template<typename ST, typename DT>
struct convert_pointer<ST const * __restrict__, DT> {
  static_assert(!std::is_pointer<ST>::value, "ST cannot be a pointer type");
  static_assert(!std::is_pointer<DT>::value, "DT cannot be a pointer type");
  using type = DT const * __restrict__;
};

template<typename ST, typename DT>
struct convert_pointer<ST *, DT> {
  static_assert(!std::is_pointer<ST>::value, "ST cannot be a pointer type");
  static_assert(!std::is_pointer<DT>::value, "DT cannot be a pointer type");
  using type = DT *;
};

template<typename ST, typename DT>
struct convert_pointer<ST * __restrict__, DT> {
  static_assert(!std::is_pointer<ST>::value, "ST cannot be a pointer type");
  static_assert(!std::is_pointer<DT>::value, "DT cannot be a pointer type");
  using type = DT * __restrict__;
};

template<typename... T>
constexpr size_t sum_sizes(std::tuple<T...> const &) {
  size_t sum = 0;
  auto a = { sum += sizeof(T)... };
  return sum;
}

namespace details {

struct apply_unary_helper_t {};

template<class T> 
T&& operator,(T&& t, apply_unary_helper_t) {
  return std::forward<T>(t); 
}

template <typename Ftor, typename Tuple, size_t... Is>
decltype(auto) apply_with_tuple_result(Ftor&& ftor, Tuple&& tuple,
std::index_sequence<Is...>) {
  auto r = std::make_tuple(
    (
      ftor(std::get<Is>(std::forward<Tuple>(tuple))),
      apply_unary_helper_t{}
    )...
  );
  return r;
}

template <typename Ftor, typename Tuple, size_t... Is>
decltype(auto) apply_with_initlist_result(Ftor&& ftor, Tuple&& tuple,
std::index_sequence<Is...>) {
  auto r = {
    (
      ftor(std::get<Is>(std::forward<Tuple>(tuple))),
      apply_unary_helper_t{}
    )...
  };
  return r;
}

} // end: namespace details

template <typename Ftor, typename Tuple>
decltype(auto) apply_with_tuple_result(Ftor&& ftor, Tuple&& tuple) {
  return details::apply_with_tuple_result(
    std::forward<Ftor>(ftor),
    std::forward<Tuple>(tuple),
    std::make_index_sequence<
      std::tuple_size<std::remove_reference_t<Tuple>>::value
    >{}
  );
}

template <typename Ftor, typename Tuple>
decltype(auto) apply_with_initlist_result(Ftor&& ftor, Tuple&& tuple) {
  return details::apply_with_initlist_result(
    std::forward<Ftor>(ftor),
    std::forward<Tuple>(tuple),
    std::make_index_sequence<
      std::tuple_size<std::remove_reference_t<Tuple>>::value
    >{}
  );
}

template <bool...> struct bool_pack;

template <bool... v>
using all_true = std::is_same<bool_pack<true, v...>, bool_pack<v..., true>>;

namespace detail {

template<typename GI, typename... T, size_t... I>
inline void assign(
  std::tuple<T * ...> const & dst, GI const dstidx,
  std::tuple<T...> const & src,
  std::index_sequence<I...>
) {
  auto r = std::tuple<T & ...>{
    (std::get<I>(dst)[dstidx] = std::get<I>(src)) ...
  };
}

template<typename GI, typename... T, size_t... I>
inline void assign(
  std::tuple<T * ...> const & dst, GI const dstidx,
  std::tuple<T const * ...> const & src, GI const srcidx,
  std::index_sequence<I...>
) {
  auto r = std::tuple<T & ...>{ (std::get<I>(dst)[dstidx] = std::get<I>(src)[srcidx]) ... };
}

template<typename GI, typename... T, size_t... I>
inline void assign(
  std::tuple<T...> & dst,
  std::tuple<T const * ...> const & src, GI const srcidx,
  std::index_sequence<I...>
) {
  auto r = std::tuple<T & ...>{
    (std::get<I>(dst) = std::get<I>(src)[srcidx]) ...
  };
}

template<typename GI, typename... T, size_t... I>
inline std::tuple<T * ...> pointers_at(
  std::tuple<T * ...> const & ptrs, GI const idx,
  std::index_sequence<I...>
) {
  return std::tuple<T * ...>(
    (&(std::get<I>(ptrs)[idx])) ...
  );
}

} // end: namespace detail

template<typename GI, typename... T>
inline void assign( // dst[dstidx] = src
  std::tuple<T * ...> const & dst, GI const dstidx,
  std::tuple<T...> const & src
) {
  #pragma forceinline
  detail::assign(dst, dstidx, src, std::make_index_sequence<sizeof...(T)>{});
}

template<typename GI, typename... T>
inline void assign( // dst[dstidx] = src[srcidx]
  std::tuple<T * ...> const & dst, GI const dstidx,
  std::tuple<T const * ...> const & src, GI const srcidx
) {
  #pragma forceinline
  detail::assign(dst, dstidx, src, srcidx, std::make_index_sequence<sizeof...(T)>{});
}

template<typename GI, typename... T>
inline void assign( // dst = src[srcidx]
  std::tuple<T...> & dst,
  std::tuple<T const * ...> const & src, GI const srcidx
) {
  #pragma forceinline
  detail::assign(dst, src, srcidx, std::make_index_sequence<sizeof...(T)>{});
}

template<typename GI, typename... T>
inline std::tuple<T * ...> pointers_at( // &ptrs[idx]
  std::tuple<T * ...> const & ptrs, GI const idx
) {
  #pragma forceinline
  return detail::pointers_at(ptrs, idx, std::make_index_sequence<sizeof...(T)>{});
}

} // end: namespace qsm
