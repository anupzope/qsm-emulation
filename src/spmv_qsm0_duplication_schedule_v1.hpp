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

#include "distribution_schedule_v1.hpp"

template<typename GI, typename DI>
class spmv_qsm0_duplication_schedule_v1;

template<typename VT>
class spmv_qsm0_duplication_buffer_v1 {
  distribution_buffer<VT> distbuffer;
  std::vector<VT *> outbufferptrs;
  
public:
  template<typename GI>
  VT * get_distribution_outbuffer(GI const part) {
    return distbuffer[part];
  }
  
  template<typename GI>
  VT const * get_distribution_outbuffer(GI const part) const {
    return distbuffer[part];
  }
  
  
  
  VT * * get_distribution_outbuffers() {
    return &outbufferptrs[0];
  }
  
  VT const * const * get_distribution_outbuffers() const {
    return &outbufferptrs[0];
  }
  
  
  
  template<typename GI>
  VT * get_distribution_tmpbuffer(GI const num, GI const part) {
    return distbuffer.tmp(num, part);
  }
  
  template<typename GI>
  VT const * get_distribution_tmpbuffer(GI const num, GI const part) const {
    return distbuffer.tmp(num, part);
  }
  
  
  
  template<typename GI, typename DI>
  void allocate(spmv_qsm0_duplication_schedule_v1<GI, DI> const & schedule) {
    GI const nparts = schedule.get_nparts();
    distbuffer.allocate(
      nparts,
      schedule.get_distribution_schedules()
    );
    outbufferptrs.resize(nparts);
    for(GI p = 0; p < nparts; ++p) {
      outbufferptrs[p] = distbuffer[p];
    }
  }
  
  
  
  void deallocate() {
    distbuffer.deallocate();
  }
};



template<typename GI, typename DI>
class spmv_qsm0_duplication_schedule_v1 {
  GI nparts;
  std::vector<distribution_schedule<GI, DI>> dschedules;
  
public:
  GI get_nparts() const {
    return nparts;
  }
  
  
  
  distribution_schedule<GI, DI> const &
  get_distribution_schedule(GI const part) const {
    return dschedules[part];
  }
  
  distribution_schedule<GI, DI> const *
  get_distribution_schedules() const {
    return &dschedules[0];
  }
  
  
  
  void generate(
    GI const nparts,
    GI const * const * copystreams, GI const * copystreamoffs,
    GI const nrows,
    typename std::make_unsigned<GI>::type const ntailbits,
    typename std::make_unsigned<GI>::type const nbitsperstage
  ) {
    this->nparts = nparts;
    dschedules.clear();
    dschedules.resize(nparts);
    
    ispace<GI> inspace(0, nrows);
    
    #pragma omp parallel for
    for(GI p = 0; p < nparts; ++p) {
      record_distribution_schedule(
        dschedules[p],
        ntailbits,
        nbitsperstage,
        1, // nparts
        &inspace,
        copystreamoffs[p+1]-copystreamoffs[p],
        copystreams[p]
      );
    }
  }
  
  
  
  template<int BufferSize, typename VT>
  void execute(
    VT const * input,
    spmv_qsm0_duplication_buffer_v1<VT> & buf
  ) {
    #pragma omp parallel for
    for(GI p = 0; p < nparts; ++p) {
      execute_distribution_schedule_using_qsm0<BufferSize>(
        buf.get_distribution_outbuffer(p),
        input,
        dschedules[p],
        buf.get_distribution_tmpbuffer(0, p),
        buf.get_distribution_tmpbuffer(1, p)
      );
    }
  }
  
  template<typename VT>
  size_t get_execution_data_volume() const {
    return distribution_schedule_helper<GI, DI>::template get_distribution_memsize<VT>(
      nparts, &dschedules[0]
    );
  }
};
