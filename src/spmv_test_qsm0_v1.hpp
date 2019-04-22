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
#include "spmv_qsm0_execution_schedule_v1.hpp"
#include "spmv_qsm0_duplication_schedule_v1.hpp"

template<typename GI, typename LI, typename DI, typename MT, typename VT>
class spmv_test_qsm0_v1 : public spmv_test<GI, MT, VT> {
protected:
  GI lookahead;
  GI lookback;
  typename std::make_unsigned<GI>::type ntailbits;
  typename std::make_unsigned<GI>::type nbitsperstage;
  
  spmv_qsm0_execution_schedule_v1<GI, LI, DI> executionschedule;
  spmv_qsm0_duplication_schedule_v1<GI, DI> copyschedule;
  std::vector<MT> qsm0a;
  spmv_qsm0_duplication_buffer_v1<VT> x0copy;
  
public:
  ~spmv_test_qsm0_v1() {
  }
  
  void set_parameters(
    GI lookahead, GI lookback,
    typename std::make_unsigned<GI>::type ntailbits,
    typename std::make_unsigned<GI>::type nbitsperstage
  ) {
    this->lookahead = lookahead;
    this->lookback = lookback;
    this->ntailbits = ntailbits;
    this->nbitsperstage = nbitsperstage;
  }
  
  virtual int preprocess() {
    printf("generating qsm0 v1 execution schedule\n");
    int retval = executionschedule.generate(
      this->nparts, &this->parts[0],
      this->matstr->nrows, this->matstr->rowoffs, this->matstr->colidxs,
      lookahead, lookback
    );
    
    if(retval) {
      this->errmsg = executionschedule.get_error_message(retval);
    } else {
      printf("generating qsm0 v1 copy schedule\n");
      copyschedule.generate(
        this->nparts,
        executionschedule.get_copystreams(),
        executionschedule.get_copystream_offs(),
        this->matstr->nrows,
        ntailbits, nbitsperstage
      );
      
      x0copy.allocate(copyschedule);
      
      printf("reordering matrix coefficient storage\n");
      qsm0a.resize(this->matstr->nnz);
      #pragma omp parallel for
      for(GI p = 0; p < this->nparts; ++p) {
        ispace<GI> const & partspace = executionschedule.get_partspace(p);
        ispace<GI> const & nnzspace = executionschedule.get_nnzspace(p);
        MT const * aptr = &this->a[nnzspace.start];
        MT * qsm0aptr = &qsm0a[nnzspace.start];
        GI const * n2o_colidx = executionschedule.get_n2o_colidx(p);
        for(GI i = 0; i < nnzspace.size; ++i) {
          qsm0aptr[i] = aptr[n2o_colidx[i]];
        }
      }
      
      //printf("generating copy stream for x0\n");
      //GI const * csoffs = schedule.get_copystream_offs();
      //GI const cssize = schedule.get_total_copystream_size();
      //GI const * csptr = schedule.get_copystream();
      //x0copy.resize(cssize);
      //#pragma omp parallel for
      //for(GI p = 0; p < nparts; ++p) {
      //  for(GI k = csoffs[p]; k < csoffs[p+1]; ++k) {
      //    x0copy[k] = x0[csptr[k]];
      //  }
      //}
    }
    
    return retval;
  }
  
  virtual void postprocess() {
    x0copy.deallocate();
  }
  
  virtual void execute_internal() {
    copyschedule.execute<10>(&this->x0[0], x0copy);
    executionschedule.execute<10>(
      this->x1, this->a, this->x0,
      x0copy.get_distribution_outbuffers()
    );
  }
  
  virtual void calculate_preliminary_metrics() {
    this->nflops = this->calculate_flop_count_per_nnz()*this->matstr->nnz;
    
    size_t dist_data_volume = copyschedule.template get_execution_data_volume<VT>();
    //size_t dist_data_volume = 
    //distribution_schedule_helper<GI, DI>::template get_distribution_memsize<VT>(
    //  this->nparts, copyschedule.get_distribution_schedules()
    //);
    
    GI total_x1size = 0;
    GI total_x0size = 0;
    GI total_asize = 0;
    GI total_copystream_size = 0;
    GI total_mainmap_size = 0;
    GI total_copymap_size = 0;
    for(GI p = 0; p < this->nparts; ++p) {
      total_x1size += executionschedule.get_partspace(p).size;
      total_x0size += executionschedule.get_partspace(p).size;
      total_asize += executionschedule.get_nnzspace(p).size;
      total_copystream_size += executionschedule.get_copystream_size(p);
      total_mainmap_size += executionschedule.get_mainmap_size(p);
      total_copymap_size += executionschedule.get_copymap_size(p);
    }
    
    size_t spmv_data_volume = 
     (total_x1size+total_x0size+total_copystream_size)*sizeof(VT) +
     total_asize*sizeof(MT) +
     (total_mainmap_size+total_copymap_size)*sizeof(LI);
    
    this->data_volume = spmv_data_volume + dist_data_volume;
  }
};
