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

#include "utility.hpp"

namespace qsm {

constexpr unsigned CLSize = 64;

constexpr unsigned CLDivMask = ~(CLSize-1);

constexpr unsigned CLDivShift = integer_log2(CLSize);

constexpr unsigned CLModMask = CLSize-1;

using clbitmask_type = uint64_t;

} // end: namespace qsm
