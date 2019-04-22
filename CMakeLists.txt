################################################################################
# Copyright 2019, Mississippi State University
# 
# This file is part of qsm-emulation.
# 
# qsm-emulation is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# qsm-emulation is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
################################################################################

cmake_minimum_required(VERSION 3.9)

project(qsm-emulation VERSION 1.0.0 LANGUAGES C CXX)

find_package(OpenMP)

#set(itt_install_dir "/usr/local/intel-2018.1/vtune_amplifier_2018.1.0.535340")
#set(itt_install_dir "/usr/local/intel-2019/intel-2019.u3/vtune_amplifier_2019.3.0.591499")
#add_library(itt STATIC IMPORTED)
#set_property(TARGET itt PROPERTY INTERFACE_INCLUDE_DIRECTORIES "${itt_install_dir}/include")
#set_property(TARGET itt PROPERTY INTERFACE_COMPILE_DEFINITIONS USING_ITT)
#set_property(TARGET itt PROPERTY IMPORTED_LOCATION "${itt_install_dir}/lib64/libittnotify.a")

add_executable(extract-matrix-structure src/mmio.c src/extract_matrix_structure.cpp)
target_include_directories(extract-matrix-structure PUBLIC src)
target_compile_features(extract-matrix-structure PUBLIC cxx_std_14)
target_compile_options(extract-matrix-structure PUBLIC -xAVX)
if (OPENMP_FOUND)
  target_link_libraries(extract-matrix-structure PUBLIC OpenMP::OpenMP_CXX)
endif()
#target_compile_options(extract-matrix-structure PUBLIC -g)

add_executable(spmv src/spmv.cpp)
target_include_directories(spmv PUBLIC src)
target_compile_features(spmv PUBLIC cxx_std_14)
target_compile_options(spmv PUBLIC -xAVX)
target_link_libraries(spmv PUBLIC rt)
if (OPENMP_FOUND)
  target_link_libraries(spmv PUBLIC OpenMP::OpenMP_CXX)
endif()
if (TARGET itt)
  target_link_libraries(spmv PUBLIC itt)
endif()
#target_compile_options(spmv PUBLIC -g)

