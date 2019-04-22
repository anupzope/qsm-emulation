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
#include "cw.hpp"

#include <vector>

template<typename GI, typename DI>
class spmv_qsm0_duplication_schedule_v2;

template<typename ... VT>
class spmv_qsm0_duplication_buffer_v2 {
  std::tuple<VT ...> * memptr;
  std::tuple<VT ...> * * * dupoutptrs;
  std::tuple<VT ...> * * distinptrs;
  distribution_buffer<VT...> distbuffer;
  std::vector<std::tuple<VT * ...>> copystreamptrs;
  std::vector<std::tuple<VT const * ...>> constcopystreamptrs;
  
public:
  spmv_qsm0_duplication_buffer_v2() : memptr{}, dupoutptrs{},
  distinptrs{} {
  }
  
  
  
  template<typename GI>
  std::tuple<VT ...> * get_duplication_outbuffer(GI const spart, GI const dpart) {
    return dupoutptrs[spart][dpart];
  }
  
  template<typename GI>
  std::tuple<VT ...> const * get_duplication_outbuffer(GI const spart, GI const dpart) const {
    return dupoutptrs[spart][dpart];
  }
  
  
  
  template<typename GI>
  std::tuple<VT ...> * get_distribution_inbuffer(GI const part) {
    return distinptrs[part];
  }
  
  template<typename GI>
  std::tuple<VT ...> * get_distribution_inbuffer(GI const part) const {
    return distinptrs[part];
  }
  
  
  
  template<typename GI>
  std::tuple<VT * ...> const & get_distribution_outbuffer(GI const part) {
    return distbuffer[part];
  }
  
  template<typename GI>
  std::tuple<VT const * ...> get_distribution_outbuffer(GI const part) const {
    return distbuffer[part];
  }
  
  
  template<typename GI>
  std::tuple<VT ...> * get_distribution_tmpbuffer(GI const num, GI const part) {
    return distbuffer.tmp(num, part);
  }
  
  template<typename GI>
  std::tuple<VT ...> const * get_distribution_tmpbuffer(GI const num, GI const part) const {
    return distbuffer.tmp(num, part);
  }
  
  
  
  template<typename GI>
  std::tuple<VT * ...> const & get_copystream(GI const part) {
    return copystreamptrs[part];
  }
  
  template<typename GI>
  std::tuple<VT const * ...> const & get_copystream(GI const part) const {
    return constcopystreamptrs[part];
  }
  
  std::tuple<VT * ...> const * get_copystreams() {
    return &copystreamptrs[0];
  }
  
  std::tuple<VT const * ...> const * get_copystreams() const {
    return &constcopystreamptrs[0];
  }
  
  
  
  template<typename GI, typename DI>
  void allocate(spmv_qsm0_duplication_schedule_v2<GI, DI> const & schedule) {
    GI const nparts = schedule.get_nparts();
    
    ispace<GI> duptspace = schedule.get_total_duplication_outspace();
    memptr = new std::tuple<VT ...>[duptspace.size];
    
    dupoutptrs = new std::tuple<VT ...>**[nparts];
    dupoutptrs[0] = new std::tuple<VT ...>*[nparts*nparts];
    for(GI p = 1; p < nparts; ++p) {
      dupoutptrs[p] = dupoutptrs[p-1] + nparts;
    }
    
    for(GI sp = 0; sp < nparts; ++sp) {
      for(GI dp = 0; dp < nparts; ++dp) {
        ispace<GI> const & spc = schedule.get_duplication_outspace(sp, dp);
        dupoutptrs[sp][dp] = &memptr[spc.start];
      }
    }
    
    distinptrs = new std::tuple<VT ...>*[nparts];
    for(GI p = 0; p < nparts; ++p) {
      ispace<GI> const & spc = schedule.get_distribution_inspace(p);
      distinptrs[p] = &memptr[spc.start];
    }
    
    distbuffer.allocate(nparts, schedule.get_distribution_schedules());
    
    copystreamptrs.resize(nparts);
    constcopystreamptrs.resize(nparts);
    for(GI p = 0; p < nparts; ++p) {
      ispace<GI> const & spc = schedule.get_copystream_space(p);
      copystreamptrs[p] = qsm::pointers_at(distbuffer[0], spc.start);
      constcopystreamptrs[p] = qsm::pointers_at(distbuffer[0], spc.start);
    }
  }
  
  void deallocate() {
    delete[] memptr;
    delete[] dupoutptrs[0];
    delete[] dupoutptrs;
    delete[] distinptrs;
    
    distbuffer.deallocate();
  }
};

template<typename VT>
class spmv_qsm0_duplication_buffer_v2<VT> {
  VT * memptr;
  VT * * * dupoutptrs;
  VT * * distinptrs;
  distribution_buffer<VT> distbuffer;
  std::vector<VT *> copystreamptrs;
  
public:
  spmv_qsm0_duplication_buffer_v2() : memptr(nullptr), dupoutptrs(nullptr),
  distinptrs(nullptr) {
  }
  
  
  
  template<typename GI>
  VT * get_duplication_outbuffer(GI const spart, GI const dpart) {
    return dupoutptrs[spart][dpart];
  }
  
  template<typename GI>
  VT const * get_duplication_outbuffer(GI const spart, GI const dpart) const {
    return dupoutptrs[spart][dpart];
  }
  
  
  
  template<typename GI>
  VT * get_distribution_inbuffer(GI const part) {
    return distinptrs[part];
  }
  
  template<typename GI>
  VT const * get_distribution_inbuffer(GI const part) const {
    return distinptrs[part];
  }
  
  
  
  template<typename GI>
  VT * get_distribution_outbuffer(GI const part) {
    return distbuffer[part];
  }
  
  template<typename GI>
  VT const * get_distribution_outbuffer(GI const part) const {
    return distbuffer[part];
  }
  
  
  
  template<typename GI>
  VT * get_distribution_tmpbuffer(GI const num, GI const part) {
    return distbuffer.tmp(num, part);
  }
  
  template<typename GI>
  VT const * get_distribution_tmpbuffer(GI const num, GI const part) const {
    return distbuffer.tmp(num, part);
  }
  
  
  
  template<typename GI>
  VT * get_copystream(GI const part) {
    return copystreamptrs[part];
  }
  
  template<typename GI>
  VT const * get_copystream(GI const part) const {
    return copystreamptrs[part];
  }
  
  VT * const * get_copystreams() {
    return &copystreamptrs[0];
  }
  
  VT const * const * get_copystreams() const {
    return &copystreamptrs[0];
  }
  
  
  
  template<typename GI, typename DI>
  void allocate(spmv_qsm0_duplication_schedule_v2<GI, DI> const & schedule) {
    GI const nparts = schedule.get_nparts();
    
    ispace<GI> duptspace = schedule.get_total_duplication_outspace();
    memptr = new VT[duptspace.size];
    
    dupoutptrs = new VT**[nparts];
    dupoutptrs[0] = new VT*[nparts*nparts];
    for(GI p = 1; p < nparts; ++p) {
      dupoutptrs[p] = dupoutptrs[p-1] + nparts;
    }
    
    for(GI sp = 0; sp < nparts; ++sp) {
      for(GI dp = 0; dp < nparts; ++dp) {
        ispace<GI> const & spc = schedule.get_duplication_outspace(sp, dp);
        dupoutptrs[sp][dp] = &memptr[spc.start];
      }
    }
    
    distinptrs = new VT*[nparts];
    for(GI p = 0; p < nparts; ++p) {
      ispace<GI> const & spc = schedule.get_distribution_inspace(p);
      distinptrs[p] = &memptr[spc.start];
    }
    
    distbuffer.allocate(nparts, schedule.get_distribution_schedules());
    
    copystreamptrs.resize(nparts);
    for(GI p = 0; p < nparts; ++p) {
      ispace<GI> const & spc = schedule.get_copystream_space(p);
      copystreamptrs[p] = &distbuffer[0][spc.start];
    }
  }
  
  void deallocate() {
    delete[] memptr;
    delete[] dupoutptrs[0];
    delete[] dupoutptrs;
    delete[] distinptrs;
    
    distbuffer.deallocate();
  }
}; // end: duplication_buffer

template<typename GI, typename DI>
class spmv_qsm0_duplication_schedule_v2 {
  GI nparts;
  std::vector<ispace<GI>> dupinspaces;
  std::vector<ispace<GI>> dupoutspaceobjs;
  std::vector<ispace<GI> *> dupoutspaces;
  std::vector<qsm::cw<4>> dupmap;
  std::vector<qsm::cw<4> *> dupmapptrs;
  std::vector<ispace<GI>> distinspaces;
  std::vector<ispace<GI>> distoutspaces;
  std::vector<distribution_schedule<GI, DI>> distschedules;
  std::vector<ispace<GI>> copystreamspaces;
  
public:
  GI get_nparts() const {
    return nparts;
  }
  
  ispace<GI> const & get_duplication_inspace(GI const part) const {
    return dupinspaces[part];
  }
  
  ispace<GI> const * get_duplication_inspaces() const {
    return &dupinspaces[0];
  }
  
  ispace<GI> get_total_duplication_inspace() const {
    return ispace<GI>(dupinspaces[0].start, dupinspaces[nparts-1].end);
  }
  
  
  
  
  ispace<GI> const * get_duplication_outspaces(GI const spart) const {
    return dupoutspaces[spart];
  }
  
  ispace<GI> const & get_duplication_outspace(GI const spart, GI const dpart) const {
    return dupoutspaces[spart][dpart];
  }
  
  ispace<GI> get_total_duplication_outspace() const {
    return ispace<GI>(dupoutspaces[0][0].start, dupoutspaces[nparts-1][nparts-1].end);
  }
  
  
  
  
  qsm::cw<4> const * get_duplication_map(GI const part) const {
    return &dupmap[dupinspaces[part].start];
  }
  
  GI get_duplication_map_size(GI const part) const {
    return dupinspaces[part].size;
  }
  
  qsm::cw<4> const * const * get_duplication_maps() const {
    return &dupmapptrs[0];
  }
  
  
  
  
  ispace<GI> const & get_distribution_inspace(GI const part) const {
    return distinspaces[part];
  }
  
  ispace<GI> const * get_distribution_inspaces() const {
    return &distinspaces[0];
  }
  
  ispace<GI> get_total_distribution_inspace() const {
    return ispace<GI>(distinspaces[0].start, distinspaces[nparts-1].end);
  }
  
  
  
  
  ispace<GI> const & get_distribution_outspace(GI const part) const {
    return distoutspaces[part];
  }
  
  ispace<GI> const * get_distribution_outspaces() const {
    return &distoutspaces[0];
  }
  
  ispace<GI> get_total_distribution_outspace() const {
    return ispace<GI>(distoutspaces[0].start, distoutspaces[nparts-1].end);
  }
  
  
  
  distribution_schedule<GI, DI> const & get_distribution_schedule(GI const part) const {
    return distschedules[part];
  }
  
  distribution_schedule<GI, DI> const * get_distribution_schedules() const {
    return &distschedules[0];
  }
  
  
  
  ispace<GI> const & get_copystream_space(GI const part) const {
    return copystreamspaces[part];
  }
  
  ispace<GI> const * get_copystream_spaces() const {
    return &copystreamspaces[0];
  }
  
  
  
  void generate(
    GI const nparts, GI const nrows,
    typename std::make_unsigned<GI>::type const ntailbits,
    typename std::make_unsigned<GI>::type const nbitsperstage,
    GI const tcssize, GI const * tcs, GI const * csoffs
  ) {
    this->nparts = nparts;
    dupinspaces.clear();
    dupoutspaceobjs.clear();
    dupoutspaces.clear();
    dupmap.clear();
    dupmapptrs.clear();
    distinspaces.clear();
    distoutspaces.clear();
    distschedules.clear();
    copystreamspaces.clear();
    
    dupinspaces.resize(nparts, ispace<GI>(0, 0));
    dupoutspaceobjs.resize(nparts*nparts, ispace<GI>(0, 0));
    dupoutspaces.resize(nparts, nullptr);
    for(GI p = 0; p < nparts; ++p) {
      dupoutspaces[p] = &dupoutspaceobjs[p*nparts];
    }
    dupmap.resize(nrows);
    for(GI r = 0; r < nrows; ++r) {
      dupmap[r].clear();
    }
    dupmapptrs.resize(nparts);
    distinspaces.resize(nparts, ispace<GI>(0, 0));
    distoutspaces.resize(nparts, ispace<GI>(0, 0));
    distschedules.resize(nparts);
    copystreamspaces.resize(nparts);
    
    for(GI p = 0; p < nparts; ++p) {
      copystreamspaces[p] = ispace<GI>(csoffs[p], csoffs[p+1]);
    }
    
    // form distribution outspaces by uniformly partitioning the index space
    // of the total copy stream
    std::vector<GI> tcspartoffs(nparts+1, 0);
    for(GI p = 0; p < nparts; ++p) {
      GI part_start, part_end;
      qsm::indexspace_partition(0, tcssize, nparts, p, &part_start, &part_end);
      tcspartoffs[p+1] = part_end;
    }
    
    // form distribution outspaces
    for(GI p = 0; p < nparts; ++p) {
      distoutspaces[p] = ispace<GI>(tcspartoffs[p], tcspartoffs[p+1]);
    }
    
    // record distribution schedules
    for(GI p = 0; p < nparts; ++p) {
      record_distribution_schedule_for_seggregated_values(
        distschedules[p],
        ntailbits, nbitsperstage,
        distoutspaces[p].size, &tcs[distoutspaces[p].start]
      );
    }
    
    // generate duplication map
    for(GI p = 0; p < nparts; ++p) {
      for(GI i = tcspartoffs[p]; i < tcspartoffs[p+1]; ++i) {
        dupmap[tcs[i]].set(p);
      }
    }
    
    // count the number of reads and writes required for each entity
    std::vector<GI> duprwcount(nrows);
    for(GI r = 0; r < nrows; ++r) {
      duprwcount[r] = dupmap[r].count() + 1;
    }
    
    // form cumulative offset of read write counts
    std::vector<GI> duprwoffs(nrows+1);
    duprwoffs[0] = 0;
    for(GI r = 0; r < nrows; ++r) {
      duprwoffs[r+1] = duprwoffs[r] + duprwcount[r];
    }
    
    // form duplication inspace offsets
    std::vector<GI> dupinspaceoffs(nparts+1);
    dupinspaceoffs[0] = 0;
    for(GI p = 0; p < nparts; ++p) {
      GI part_start, part_end;
      qsm::indexspace_partition(0, duprwoffs[nrows], nparts, p, &part_start, &part_end);
      dupinspaceoffs[p+1] = part_end;
    }
    for(GI r = 0, p = 1; r < nrows; ++r) {
      if(dupinspaceoffs[p] <= duprwoffs[r+1]) {
        dupinspaceoffs[p++] = r+1;
      }
    }
    
    // form duplication inspaces
    for(GI sp = 0; sp < nparts; ++sp) {
      dupinspaces[sp] = ispace<GI>(dupinspaceoffs[sp], dupinspaceoffs[sp+1]);
      dupmapptrs[sp] = &dupmap[dupinspaceoffs[sp]];
    }
    
    // form source to destination duplication count matrix
    std::vector<std::vector<GI>> dupoutmat(nparts, std::vector<GI>(nparts, 0));
    for(GI p = 0; p < nparts; ++p) {
      for(GI i = dupinspaces[p].start; i < dupinspaces[p].end; ++i) {
        qsm::cw<4>::iterator m = dupmap[i].begin();
        while(m.has_next()) {
          int dp = m.next();
          dupoutmat[p][dp]++;
        }
      }
    }
    
    // form duplication outspaces
    // form distribution outspaces
    for(GI dp = 0, off = 0; dp < nparts; ++dp) {
      for(GI sp = 0; sp < nparts; ++sp) {
        dupoutspaces[sp][dp] = ispace<GI>(off, off + dupoutmat[sp][dp]);
        off += dupoutmat[sp][dp];
      }
      distinspaces[dp] = ispace<GI>(
        dupoutspaces[0][dp].start, dupoutspaces[nparts-1][dp].end
      );
    }
  }
  
  
  
  template<int BufferSize, typename ... VT>
  void execute(
    std::tuple<VT const * ...> const & in,
    spmv_qsm0_duplication_buffer_v2<VT...> & out
  ) {
//#define SPMV_PACK_QSM0_DUPLICATION_SCHEDULE_V2_EXE_TEMPORAL
#ifdef SPMV_PACK_QSM0_DUPLICATION_SCHEDULE_V2_EXE_TEMPORAL
    #pragma omp parallel for
    for(GI sp = 0; sp < nparts; ++sp) {
      ispace<GI> const & inspc = get_duplication_inspace(sp);
      std::tuple<VT const * ...> inptr = qsm::pointers_at(in, inspc.start);
      qsm::cw<4> const * map = get_duplication_map(sp);
      std::tuple<VT ...> * outptrs[0x1u << 4];
      GI outoffs[0x1u << 4];
      for(GI dp = 0; dp < nparts; ++dp) {
        outptrs[dp] = out.get_duplication_outbuffer(sp, dp);
        outoffs[dp] = 0;
      }
      for(GI i = 0; i < inspc.size; ++i) {
        qsm::cw<4>::iterator it = map[i].begin();
        while(it.has_next()) {
          int p = it.next();
          //*(outptrs[p]) = inptr[i];
          qsm::assign(outptrs[p][outoffs[p]], inptr, i);
          outoffs[p]++;
        }
      }
    }
    
    #pragma omp parallel for
    for(GI p = 0; p < nparts; ++p) {
      execute_distribution_schedule_in_pack(
        out.get_distribution_outbuffer(p),
        out.get_distribution_inbuffer(p),
        get_distribution_schedule(p),
        out.get_distribution_tmpbuffer(0, p),
        out.get_distribution_tmpbuffer(1, p)
      );
    }
#else
    #pragma omp parallel for
    for(GI sp = 0; sp < nparts; ++sp) {
      ispace<GI> const & inspc = get_duplication_inspace(sp);
      std::tuple<VT const * ...> inptr = qsm::pointers_at(in, inspc.start);
      qsm::cw<4> const * map = get_duplication_map(sp);
      
      qsm::cgwqueue<GI, std::tuple<VT ...>, BufferSize> wq[0x1u << 4];
      std::tuple<VT ...> * bufptrs[0x1u << 4];
      int bufidx[0x1u << 4];
      
      for(GI dp = 0; dp < nparts; ++dp) {
        wq[dp].configure(out.get_duplication_outbuffer(sp, dp));
      }
      
      for(GI i = 0; i < inspc.size; i+=BufferSize) {
        for(GI dp = 0; dp < nparts; ++dp) {
          bufptrs[dp] = wq[dp].get();
          bufidx[dp] = 0;
        }
        
        int niter = std::min(inspc.size-i, BufferSize);
        #pragma loop_count max(BufferSize)
        for(int k = 0; k < niter; ++k) {
          qsm::cw<4>::iterator it = map[i+k].begin();
          while(it.has_next()) {
            int p = it.next();
            //bufptrs[p][bufidx[p]] = inptr[i+k];
            qsm::assign(bufptrs[p][bufidx[p]], inptr, i+k);
            bufidx[p]++;
          }
        }
        
        for(GI dp = 0; dp < nparts; ++dp) {
          wq[dp].flush(bufidx[dp]);
        }
      }
      
      for(GI dp = 0; dp < nparts; ++dp) {
        wq[dp].deactivate();
      }
    }
    
    #pragma omp parallel for
    for(GI p = 0; p < nparts; ++p) {
      execute_distribution_schedule_in_pack_using_qsm0<BufferSize>(
        out.get_distribution_outbuffer(p),
        out.get_distribution_inbuffer(p),
        get_distribution_schedule(p),
        out.get_distribution_tmpbuffer(0, p),
        out.get_distribution_tmpbuffer(1, p)
      );
    }
#endif // #ifdef SPMV_PACK_QSM0_DUPLICATION_SCHEDULE_V2_EXE_TEMPORAL
  }
  
  template<int BufferSize, typename VT>
  void execute(VT const * in, spmv_qsm0_duplication_buffer_v2<VT> & out) {
//#define SPMV_QSM0_DUPLICATION_SCHEDULE_V2_EXE_TEMPORAL
#ifdef SPMV_QSM0_DUPLICATION_SCHEDULE_V2_EXE_TEMPORAL
    // Temporal version for duplication
    #pragma omp parallel for
    for(GI sp = 0; sp < nparts; ++sp) {
      ispace<GI> const & inspc = get_duplication_inspace(sp);
      VT const * inptr = &in[inspc.start];
      qsm::cw<4> const * map = get_duplication_map(sp);
      VT * outptrs[0x1u << 4];
      for(GI dp = 0; dp < nparts; ++dp) {
        outptrs[dp] = out.get_duplication_outbuffer(sp, dp);
      }
      for(GI i = 0; i < inspc.size; ++i) {
        qsm::cw<4>::iterator it = map[i].begin();
        while(it.has_next()) {
          int p = it.next();
          *(outptrs[p]) = inptr[i];
          outptrs[p]++;
        }
      }
    }
    
    #pragma omp parallel for
    for(GI p = 0; p < nparts; ++p) {
      execute_distribution_schedule(
        out.get_distribution_outbuffer(p),
        out.get_distribution_inbuffer(p),
        get_distribution_schedule(p),
        out.get_distribution_tmpbuffer(0, p),
        out.get_distribution_tmpbuffer(1, p)
      );
    }
#else
    #pragma omp parallel for
    for(GI sp = 0; sp < nparts; ++sp) {
      ispace<GI> const & inspc = get_duplication_inspace(sp);
      VT const * inptr = &in[inspc.start];
      qsm::cw<4> const * map = get_duplication_map(sp);
      
      qsm::cgwqueue<GI, VT, BufferSize> wq[0x1u << 4];
      VT * bufptrs[0x1u << 4];
      int bufidx[0x1u << 4];
      
      for(GI dp = 0; dp < nparts; ++dp) {
        wq[dp].configure(out.get_duplication_outbuffer(sp, dp));
      }
      
      for(GI i = 0; i < inspc.size; i+=BufferSize) {
        for(GI dp = 0; dp < nparts; ++dp) {
          bufptrs[dp] = wq[dp].get();
          bufidx[dp] = 0;
        }
        
        int niter = std::min(inspc.size-i, BufferSize);
        #pragma loop_count max(BufferSize)
        for(int k = 0; k < niter; ++k) {
          qsm::cw<4>::iterator it = map[i+k].begin();
          while(it.has_next()) {
            int p = it.next();
            bufptrs[p][bufidx[p]] = inptr[i+k];
            bufidx[p]++;
          }
        }
        
        for(GI dp = 0; dp < nparts; ++dp) {
          wq[dp].flush(bufidx[dp]);
        }
      }
      
      for(GI dp = 0; dp < nparts; ++dp) {
        wq[dp].deactivate();
      }
    }
    
    #pragma omp parallel for
    for(GI p = 0; p < nparts; ++p) {
      execute_distribution_schedule_using_qsm0<BufferSize>(
        out.get_distribution_outbuffer(p),
        out.get_distribution_inbuffer(p),
        get_distribution_schedule(p),
        out.get_distribution_tmpbuffer(0, p),
        out.get_distribution_tmpbuffer(1, p)
      );
    }
#endif // #ifdef SPMV_QSM0_DUPLICATION_SCHEDULE_V2_EXE_TEMPORAL
  }
  
  template<typename VT>
  size_t get_execution_data_volume() const {
    size_t dupinsize = 0, dupoutsize = 0, dupmapsize = 0;
    for(GI sp = 0; sp < nparts; ++sp) {
      dupinsize += get_duplication_inspace(sp).size;
      for(GI dp = 0; dp < nparts; ++dp) {
        dupoutsize += get_duplication_outspace(sp, dp).size;
      }
      dupmapsize += get_duplication_map_size(sp);
    }
    
    size_t vol1 = (dupinsize+dupoutsize)*sizeof(VT)+dupmapsize*sizeof(qsm::cw<4>);
    
    size_t distinsize, distoutsize, distmapsize;
    size_t vol2 = distribution_schedule_helper<GI, DI>::template get_distribution_memsize<VT>(
      nparts, &distschedules[0], &distinsize, &distoutsize, &distmapsize
    );
    
    //printf("insize = %lu, outsize = %lu, mapsize = %lu\n",
    //  insize/1024.0/1024.0, outsize/1024.0/1024.0, mapsize/1024.0/1024.0);
    //printf("duplication volume = %lf, distribution volume = %lf\n",
    //  vol1/1024.0/1024.0, vol2/1024.0/1024.0);
    
    return vol1 + vol2;
  }
  
  template<typename ... VT>
  size_t get_execution_data_volume() const {
    size_t sz1 = qsm::sum_sizes(std::tuple<VT ...>());
    size_t sz2 = sizeof(std::tuple<VT ...>);
    
    size_t dupinsize = 0, dupoutsize = 0, dupmapsize = 0;
    for(GI sp = 0; sp < nparts; ++sp) {
      dupinsize += get_duplication_inspace(sp).size;
      for(GI dp = 0; dp < nparts; ++dp) {
        dupoutsize += get_duplication_outspace(sp, dp).size;
      }
      dupmapsize += get_duplication_map_size(sp);
    }
    
    size_t vol1 = dupinsize*sz1 + dupoutsize*sz2 + dupmapsize*sizeof(qsm::cw<4>);
    
    size_t distinsize, distoutsize, distmapsize;
    size_t vol2 = distribution_schedule_helper<GI, DI>::template get_distribution_memsize3<VT...>(
      nparts, &distschedules[0], &distinsize, &distoutsize, &distmapsize
    );
    
    //printf("insize = %lu, outsize = %lu, mapsize = %lu\n",
    //  insize/1024.0/1024.0, outsize/1024.0/1024.0, mapsize/1024.0/1024.0);
    //printf("duplication volume = %lf, distribution volume = %lf\n",
    //  vol1/1024.0/1024.0, vol2/1024.0/1024.0);
    
    return vol1 + vol2;
  }
}; // end: spmv_qsm0_duplication_schedule_v2
