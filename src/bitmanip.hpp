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

template<typename T>
T bit1atn(int index) {
  static_assert(false, "Not implemented");
}

template<typename T>
inline T setbit(T v, int index) {
  static_assert(false, "Not implemented");
}

template<typename T>
inline T setzero() {
  static_assert(false, "Not implemented");
}

template<typename T>
inline bool testzero(T v) {
  static_assert(false, "Not implemented");
}

template<typename T>
inline bool testnotzero(T v) {
  static_assert(false, "Not implemented");
}

template<typename T>
inline bool testnotequal(T a, T b) {
  static_assert(false, "Not implemented");
}

template<typename T>
inline bool testequal(T a, T b) {
  static_assert(false, "Not implemented");
}

template<typename T>
inline int bit1count(T v) {
  static_assert(false, "Not implemented");
}

template<typename T>
inline int ls1bindex(T v) {
  static_assert(false, "Not implemented");
}

template<typename T>
inline T log2(T v) {
  static_assert(false, "Not implemented");
}

} // end: namespace qsm

#if defined(__AVX__)
#include "bitmanip_avx.hpp"
#else
#error "Could not find platform specific implementation of bitmanip.hpp"
#endif

