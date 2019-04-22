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

#include "write_queue.hpp"

template<typename GI, typename LI, typename DI>
class spmv_qsm0_execution_schedule_v1 {
  GI nparts;
  std::vector<ispace<GI>> partspaces;
  std::vector<ispace<GI>> nnzspaces;
  
  std::vector<GI> copystreams;
  std::vector<GI> copystreamoffs;
  std::vector<GI *> copystreamptrs;
  
  std::vector<LI> mainmaps;
  std::vector<GI> mainmapoffs;
  std::vector<LI *> mainmapptrs;
  
  std::vector<LI> copymaps;
  std::vector<GI> copymapoffs;
  std::vector<LI *> copymapptrs;
  
  std::vector<GI> permcolidxs;
  std::vector<GI> permcolidxoffs;
  std::vector<GI *> permcolidxptrs;
  
public:
  void reset(GI const nparts) {
    this->nparts = nparts;
    partspaces.clear();
    partspaces.resize(nparts);
    nnzspaces.clear();
    nnzspaces.resize(nparts);
    
    copystreams.clear();
    copystreamoffs.clear();
    copystreamoffs.resize(nparts+1, 0);
    copystreamptrs.clear();
    copystreamptrs.resize(nparts, nullptr);
    
    mainmaps.clear();
    mainmapoffs.clear();
    mainmapoffs.resize(nparts+1, 0);
    mainmapptrs.clear();
    mainmapptrs.resize(nparts, nullptr);
    
    copymaps.clear();
    copymapoffs.clear();
    copymapoffs.resize(nparts+1, 0);
    copymapptrs.clear();
    copymapptrs.resize(nparts, nullptr);
    
    permcolidxs.clear();
    permcolidxoffs.clear();
    permcolidxoffs.resize(nparts+1, 0);
    permcolidxptrs.clear();
    permcolidxptrs.resize(nparts, nullptr);
  }
  
  void set_partspace(GI const part, ispace<GI> const & space) {
    partspaces[part] = space;
  }
  
  void set_nnzspace(GI const part, ispace<GI> const & space) {
    nnzspaces[part] = space;
  }
  
  void set_copystream(GI const part, GI const size, GI const * values) {
    GI oldsize = copystreamoffs[part+1]-copystreamoffs[part];
    if(oldsize > 0) {
      auto startit = copystreams.begin()+copystreamoffs[part];
      auto endit = copystreams.begin()+copystreamoffs[part+1];
      copystreams.erase(startit, endit);
    }
    
    auto startit = copystreams.begin()+copystreamoffs[part];
    copystreams.insert(startit, values, values+size);
    
    for(GI p = part; p < nparts; ++p) {
      copystreamoffs[p+1] -= oldsize;
      copystreamoffs[p+1] += size;
    }
    
    for(GI p = 0; p < nparts; ++p) {
      copystreamptrs[p] = &copystreams[copystreamoffs[p]];
    }
  }
  
  void set_mainmap(GI const part, GI const size, LI const * values) {
    GI oldsize = mainmapoffs[part+1]-mainmapoffs[part];
    if(oldsize > 0) {
      auto startit = mainmaps.begin()+mainmapoffs[part];
      auto endit = mainmaps.begin()+mainmapoffs[part+1];
      mainmaps.erase(startit, endit);
    }
    
    auto startit = mainmaps.begin()+mainmapoffs[part];
    mainmaps.insert(startit, values, values+size);
    
    for(GI p = part; p < nparts; ++p) {
      mainmapoffs[p+1] -= oldsize;
      mainmapoffs[p+1] += size;
    }
    
    for(GI p = 0; p < nparts; ++p) {
      mainmapptrs[p] = &mainmaps[mainmapoffs[p]];
    }
  }
  
  void set_copymap(GI const part, GI const size, LI const * values) {
    GI oldsize = copymapoffs[part+1]-copymapoffs[part];
    if(oldsize > 0) {
      auto startit = copymaps.begin()+copymapoffs[part];
      auto endit = copymaps.begin()+copymapoffs[part+1];
      copymaps.erase(startit, endit);
    }
    
    auto startit = copymaps.begin()+copymapoffs[part];
    copymaps.insert(startit, values, values+size);
    
    for(GI p = part; p < nparts; ++p) {
      copymapoffs[p+1] -= oldsize;
      copymapoffs[p+1] += size;
    }
    
    for(GI p = 0; p < nparts; ++p) {
      copymapptrs[p] = &copymaps[copymapoffs[p]];
    }
  }
  
  void set_n2o_colidx(GI const part, GI const size, GI const * values) {
    GI oldsize = permcolidxoffs[part+1]-permcolidxoffs[part];
    if(oldsize > 0) {
      auto startit = permcolidxs.begin()+permcolidxoffs[part];
      auto endit = permcolidxs.begin()+permcolidxoffs[part+1];
      permcolidxs.erase(startit, endit);
    }
    
    auto startit = permcolidxs.begin()+permcolidxoffs[part];
    permcolidxs.insert(startit, values, values+size);
    
    for(GI p = part; p < nparts; ++p) {
      permcolidxoffs[p+1] -= oldsize;
      permcolidxoffs[p+1] += size;
    }
    
    for(GI p = 0; p < nparts; ++p) {
      permcolidxptrs[p] = &permcolidxs[permcolidxoffs[p]];
    }
  }
  
public:
  GI get_nparts() const {
    return nparts;
  }
  
  
  
  ispace<GI> const & get_partspace(GI const part) const {
    return partspaces[part];
  }
  
  ispace<GI> const * get_partspaces() const {
    return &partspaces[0];
  }
  
  
  
  ispace<GI> const & get_nnzspace(GI const part) const {
    return nnzspaces[part];
  }
  
  ispace<GI> const * get_nnzspaces() const {
    return &nnzspaces[0];
  }
  
  
  
  GI const * get_copystream(GI const part) const {
    return &copystreams[copystreamoffs[part]];
  }
  
  GI get_copystream_size(GI const part) const {
    return copystreamoffs[part+1]-copystreamoffs[part];
  }
  
  GI const * const * get_copystreams() const {
    return &copystreamptrs[0];
  }
  
  GI const * get_total_copystream() const {
    return &copystreams[0];
  }
  
  GI get_total_copystream_size() const {
    return copystreams.size();
  }
  
  GI const * get_copystream_offs() const {
    return &copystreamoffs[0];
  }
  
  
  
  LI const * get_mainmap(GI const part) const {
    return &mainmaps[mainmapoffs[part]];
  }
  
  GI get_mainmap_size(GI const part) const {
    return mainmapoffs[part+1]-mainmapoffs[part];
  }
  
  LI const * const * get_mainmaps() const {
    return &mainmapptrs[0];
  }
  
  LI const * get_total_mainmap() const {
    return &mainmaps[0];
  }
  
  GI get_total_mainmap_size() const {
    return mainmaps.size();
  }
  
  GI const * get_mainmap_offs() const {
    return &mainmapoffs[0];
  }
  
  
  
  LI const * get_copymap(GI const part) const {
    return &copymaps[copymapoffs[part]];
  }
  
  GI get_copymap_size(GI const part) const {
    return copymapoffs[part+1]-copymapoffs[part];
  }
  
  LI const * const * get_copymaps() const {
    return &copymapptrs[0];
  }
  
  LI const * get_total_copymap() const {
    return &copymap[0];
  }
  
  GI get_total_copymap_size() const {
    return copymaps.size();
  }
  
  GI const * get_copymap_offs() const {
    return &copymapoffs[0];
  }
  
  
  
  GI const * get_n2o_colidx(GI const part) const {
    return &permcolidxs[permcolidxoffs[part]];
  }
  
  GI get_n2o_colidx_size(GI const part) const {
    return permcolidxoffs[part+1]-permcolidxoffs[part];
  }
  
  GI const * const * get_n2o_colidxs() const {
    return &permcolidxptrs[0];
  }
  
  GI const * get_total_n2o_colidx() const {
    return &permcolidxs[0];
  }
  
  GI get_total_n2o_colidx_size() const {
    return permcolidxs.size();
  }
  
  GI const * get_n2o_colidx_offs() const {
    return &permcolidxoffs[0];
  }
  
  char const * get_error_message(int errcode) const {
    static char const * msgs[] = {
      "success",
      "node degree too large for LI type"
    };
    
    return msgs[errcode];
  }
  
  int generate(
    GI const nparts, ispace<GI> const * parts,
    GI const nrows, GI const * rowoffs, GI const * colidxs,
    GI const advance, GI const lookback
  ) {
    int err = 0;
    
    reset(nparts);
    
    #pragma omp parallel for
    for(GI t = 0; t < nparts; ++t) {
      int local_err = 0;
      
      #pragma omp atomic read
      local_err = err;
      
      if(local_err) {
        continue;
      }
      
      GI const part_low = parts[t].start;
      GI const part_high = parts[t].end;
      GI const nnzstart = rowoffs[part_low];
      GI const nnzend = rowoffs[part_high];
      GI const localnnz = nnzend-nnzstart;
      
      std::vector<GI> copystream;
      std::vector<LI> mainmap;
      std::vector<LI> copymap;
      std::vector<GI> colidxn2o(localnnz);
      
      for(GI r = part_low; r < part_high; ++r) {
        // Determine extents of queue window
        GI high = r + advance;
        if(high > part_high) high = part_high;
        GI low = part_low;
        if(high > (part_low+lookback))
          low = high - lookback;
        
        GI nmainconn = 0, ncopylbconn = 0, ncopysqconn = 0;
        // 0: main, 1: copy lookback, 2: copy sequential
        // original offset in colidxs
        // original column idx
        // look-back index
        std::vector<std::tuple<GI, GI, GI, LI>> colidxinfo(rowoffs[r+1]-rowoffs[r]);
        for(GI j = rowoffs[r], k = 0; j < rowoffs[r+1]; ++j, ++k) {
          GI c = colidxs[j];
          if(c >= low && c < high) {
            colidxinfo[k] = std::tuple<GI, GI, GI, LI>(0, j, c, c-r);
            ++nmainconn;
          } else {
            auto cit = copystream.rbegin();
            bool found = false;
            GI dist = 0;
            while(cit != copystream.rend() && dist < lookback) {
              ++dist;
              if(*cit == c) { found = true; break; }
              ++cit;
            }
            if(found) {
              colidxinfo[k] = std::tuple<GI, GI, GI, LI>(1, j, c, dist);
              ++ncopylbconn;
            } else {
              colidxinfo[k] = std::tuple<GI, GI, GI, LI>(2, j, c, 0);
              ++ncopysqconn;
            }
          }
        }
        
        if(nmainconn < 0 || nmainconn > std::numeric_limits<LI>::max()) {
          local_err = 1;
          break;
        }
        if(ncopylbconn < 0 || ncopylbconn > std::numeric_limits<LI>::max()) {
          local_err = 1;
          break;
        }
        if(ncopysqconn < 0 || ncopysqconn > std::numeric_limits<LI>::max()) {
          local_err = 1;
          break;
        }
        
        std::sort(colidxinfo.begin(), colidxinfo.end());
        
        mainmap.push_back(nmainconn);
        copymap.push_back(ncopysqconn);
        copymap.push_back(ncopylbconn);
        
        for(GI j = rowoffs[r], k = 0; j < rowoffs[r+1]; ++j, ++k) {
          GI type = std::get<0>(colidxinfo[k]);
          GI oj = std::get<1>(colidxinfo[k]);
          GI c = std::get<2>(colidxinfo[k]);
          LI idx = std::get<3>(colidxinfo[k]);
          switch(type) {
          case 0:
            mainmap.push_back(idx);
            break;
          case 1:
            copymap.push_back(idx);
            break;
          case 2:
            copystream.push_back(c);
            break;
          }
          colidxn2o[j-nnzstart] = oj-nnzstart;
        }
      } // end: r loop
      
      if(local_err) {
        #pragma omp atomic write
        err = local_err;
      } else {
        #pragma omp critical
        {
          set_partspace(t, ispace<GI>(part_low, part_high));
          set_nnzspace(t, ispace<GI>(nnzstart, nnzend));
          set_copystream(t, copystream.size(), &copystream[0]);
          set_mainmap(t, mainmap.size(), &mainmap[0]);
          set_copymap(t, copymap.size(), &copymap[0]);
          set_n2o_colidx(t, localnnz, &colidxn2o[0]);
        }
      }
    } // end: p loop
    
    if(err) reset(0);
    
    return err;
  }
  
  template<int BufferSize, typename MT, typename ... VT>
  void execute(
    std::tuple<VT * ...> const & x1,
    MT const * a,
    std::tuple<VT const * ...> const &x0,
    std::tuple<VT const * ...> const * x0copy
  ) {
    //constexpr int BufferSize = 10;
    
    GI const nparts = this->get_nparts();
    GI const * csoffs = this->get_copystream_offs();
    
    #pragma omp parallel for
    for(GI p = 0; p < nparts; ++p) {
      ispace<GI> const & partspace = this->get_partspace(p);
      ispace<GI> const & nnzspace = this->get_nnzspace(p);
      
      std::tuple<VT * ...> x1ptr = qsm::pointers_at(x1, partspace.start);
      MT const * aptr = &a[nnzspace.start];
      std::tuple<VT const * ...> x0ptr = qsm::pointers_at(x0, partspace.start);
      std::tuple<VT const * ...> x0copyptr = x0copy[p];
      LI const * mmap = this->get_mainmap(p);
      LI const * cmap = this->get_copymap(p);
      
      qsm::cgwqueue_pack<BufferSize, GI, VT...> x1q;
      x1q.configure(x1ptr);
      std::tuple<VT * ...> bufptr;
      
      for(GI r = 0; r < partspace.size; r+=BufferSize) {
        bufptr = x1q.get();
        
        int niter = std::min(partspace.size-r, BufferSize);
        
        #pragma loop_count max(BufferSize)
        for(int k = 0; k < niter; ++k) {
          std::tuple<VT ...> sum;
          spmv_pack_helper::assign(sum, 0.0);
          
          LI nmainconn = *mmap++;
          for(GI i = 0; i < nmainconn; ++i) {
            LI idx = *mmap++;
            spmv_pack_helper::multiply_add(sum, *aptr++, x0ptr, r+k+idx);
          }
          
          LI ncopysqconn = *cmap++;
          LI ncopylbconn = *cmap++;
          for(GI i = 0; i < ncopylbconn; ++i) {
            LI idx = *cmap++;
            spmv_pack_helper::multiply_add(sum, *aptr++, x0copyptr, -idx);
          }
          
          for(GI i = 0; i < ncopysqconn; ++i) {
            spmv_pack_helper::multiply_add(sum, *aptr++, x0copyptr, i);
          }
          
          x0copyptr = qsm::pointers_at(x0copyptr, ncopysqconn);
          
          qsm::assign(bufptr, k, sum);
        } // end: k loop
        
        x1q.flush(niter);
      } // end: r loop
      
      x1q.deactivate();
    } // end: p loop
  }
  
  template<int BufferSize, typename MT, typename VT>
  void execute(
    VT * x1, MT const * a, VT const * x0,
    VT const * const * x0copy
  ) {
    //constexpr int BufferSize = 10;
    
    GI const nparts = this->get_nparts();
    GI const * csoffs = this->get_copystream_offs();
    
    #pragma omp parallel for
    for(GI p = 0; p < nparts; ++p) {
      ispace<GI> const & partspace = this->get_partspace(p);
      ispace<GI> const & nnzspace = this->get_nnzspace(p);
      
      VT * x1ptr = &x1[partspace.start];
      MT const * aptr = &a[nnzspace.start];
      VT const * x0ptr = &x0[partspace.start];
      VT const * x0copyptr = x0copy[p];
      LI const * mmap = this->get_mainmap(p);
      LI const * cmap = this->get_copymap(p);
      
      qsm::cgwqueue<GI, VT, BufferSize> x1q;
      x1q.configure(x1ptr);
      
      for(GI r = 0; r < partspace.size; r+=BufferSize) {
        auto bufptr = x1q.get();
        
        int niter = std::min(partspace.size-r, BufferSize);
        
        #pragma loop_count max(BufferSize)
        for(int k = 0; k < niter; ++k) {
          VT sum = 0.0;
          
          LI nmainconn = *mmap++;
          for(GI i = 0; i < nmainconn; ++i) {
            LI idx = *mmap++;
            sum += (*aptr++) * x0ptr[r+k+idx];
          }
          
          LI ncopysqconn = *cmap++;
          LI ncopylbconn = *cmap++;
          for(GI i = 0; i < ncopylbconn; ++i) {
            LI idx = *cmap++;
            sum += (*aptr++) * x0copyptr[-idx];
          }
          
          for(GI i = 0; i < ncopysqconn; ++i) {
            sum += (*aptr++) * (*x0copyptr++);
          }
          
          bufptr[k] = sum;
        } // end: k loop
        
        x1q.flush(niter);
      } // end: r loop
      
      x1q.deactivate();
    } // end: p loop
  }
}; // end: spmv_qsm0_execution_schedule_v1
