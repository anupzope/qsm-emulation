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

#include <cstdint>

namespace qsm {

constexpr uint8_t debruijn8_2n = (uint8_t)0x1Du;

constexpr uint8_t debruijn8table_2n[8] = {
  0, 1, 6, 2, 7, 5, 4, 3
};

constexpr uint16_t debruijn16_2n = (uint16_t)0x0D2Fu;

constexpr uint16_t debruijn16table_2n[16] = {
  0, 1, 8, 2, 6, 9, 3, 11, 15, 7, 5, 10, 14, 4, 13, 12
};

constexpr uint32_t debruijn32_2n = 0x077CB531u;

constexpr uint32_t debruijn32table_2n[32] = {
  0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
  31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
};

constexpr uint32_t debruijn32_2nm1 = 0x07C4ACDDu;

constexpr uint32_t debruijn32table_2nm1[32] = {
  0, 9, 1, 10, 13, 21, 2, 29, 11, 14, 16, 18, 22, 25, 3, 30,
  8, 12, 20, 28, 15, 17, 24, 7, 19, 27, 23, 6, 26, 5, 4, 31
};

constexpr uint64_t debruijn64_2n = 0x022FDD63CC95386Dul;

constexpr uint64_t debruijn64table_2n[64] = {
  0,  1,  2, 53,  3,  7, 54, 27,
  4, 38, 41,  8, 34, 55, 48, 28,
  62,  5, 39, 46, 44, 42, 22,  9,
  24, 35, 59, 56, 49, 18, 29, 11,
  63, 52,  6, 26, 37, 40, 33, 47,
  61, 45, 43, 21, 23, 58, 17, 10,
  51, 25, 36, 32, 60, 20, 57, 16,
  50, 31, 19, 15, 30, 14, 13, 12,
};

constexpr uint64_t debruijn64_2nm1 = 0x03F79D71B4CB0A89ul;

constexpr uint64_t debruijn64table_2nm1[64] = {
  0 , 47, 1 , 56, 48, 27, 2 , 60,
  57, 49, 41, 37, 28, 16, 3 , 61,
  54, 58, 35, 52, 50, 42, 21, 44,
  38, 32, 29, 23, 17, 11, 4 , 62,
  46, 55, 26, 59, 40, 36, 15, 53,
  34, 51, 20, 43, 31, 22, 10, 45,
  25, 39, 14, 33, 19, 30, 9 , 24,
  13, 18, 8 , 12, 7 , 6 , 5 , 63,
};

} // end: namespace qsm
