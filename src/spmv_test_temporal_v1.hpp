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

#include "spmv_test.hpp"

template<typename GI, typename MT, typename VT>
class spmv_test_temporal_v1 : public spmv_test<GI, MT, VT> {
public:
  ~spmv_test_temporal_v1() {
  }
  
  virtual void execute_internal() {
    #pragma omp parallel for
    for(GI p = 0; p < this->nparts; ++p) {
      for(GI r = this->parts[p].start; r < this->parts[p].end; ++r) {
        VT sum = 0;
        for(GI i = this->matstr->rowoffs[r]; i < this->matstr->rowoffs[r+1]; ++i) {
          sum += this->a[i]*this->x0[this->matstr->colidxs[i]];
        }
        this->x1[r] = sum;
      }
    }
  }
  
  virtual void calculate_preliminary_metrics() {
    this->nflops = spmv_traits<MT, VT>::get_flop_count_per_nnz()*this->matstr->nnz;
    
    this->data_volume = this->matstr->nnz*sizeof(MT) +
      this->matstr->nrows*sizeof(VT) +
      this->matstr->ncols*sizeof(VT) +
      (this->matstr->nrows+1)*sizeof(GI) + this->matstr->nnz*sizeof(GI);
  }
};

template<typename GI, typename MT, typename ... VT>
class spmv_pack_test_temporal_v1 : public spmv_pack_test<GI, MT, VT ...> {
public:
  ~spmv_pack_test_temporal_v1() {
  }
  
  virtual void execute_internal() {
    #pragma omp parallel for
    for(GI p = 0; p < this->nparts; ++p) {
      for(GI r = this->parts[p].start; r < this->parts[p].end; ++r) {
        std::tuple<VT ...> sum;
        spmv_pack_helper::assign(sum, 0.0);
        for(GI i = this->matstr->rowoffs[r]; i < this->matstr->rowoffs[r+1]; ++i) {
          spmv_pack_helper::multiply_add(sum, this->a[i], this->x0, this->matstr->colidxs[i]);
        }
        spmv_pack_helper::assign(this->x1, r, sum);
      }
    }
  }
  
  virtual void calculate_preliminary_metrics() {
    this->nflops = calculate_flop_count_per_nnz(
      std::make_index_sequence<sizeof...(VT)>{})*this->matstr->nnz;
    
    size_t sz = qsm::sum_sizes(std::tuple<VT...>());
    this->data_volume = this->matstr->nnz*sizeof(MT) +
      this->matstr->nrows*sz +
      this->matstr->ncols*sz +
      (this->matstr->nrows+1)*sizeof(GI) + this->matstr->nnz*sizeof(GI);
  }
};
