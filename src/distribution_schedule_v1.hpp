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

#include <utility.hpp>
#include <read_queue.hpp>
#include <write_queue.hpp>

#include <algorithm>
#include <set>
#include <vector>

template<typename GI, typename DI>
class distribution_schedule {
private:
  using UGI = typename std::make_unsigned<GI>::type;
  
public:
  using index_type = GI;
  using uindex_type = UGI;
  using map_type = DI;
  
private:
  UGI npushstages;
  UGI ntotalbits;
  UGI ntailbits;
  UGI nbitsperstage;
  UGI nbuckets;
  std::vector<ispace<GI>> input_parts;
  std::vector<std::vector<ispace<UGI>>> push_parts;
  std::vector<std::vector<ispace<UGI>>> push_subparts;
  std::vector<std::vector<DI>> push_maps;
  std::vector<std::vector<UGI>> push_mapoffs;
  std::vector<ispace<UGI>> pull_parts;
  std::vector<UGI> pull_maps;
  std::vector<UGI> pull_mapoffs;
  
public:
  distribution_schedule()
  : npushstages(0), ntotalbits(0), ntailbits(0), nbitsperstage(0), nbuckets(0) {
  }
  
  distribution_schedule(
    UGI npushstages, UGI ntotalbits, UGI ntailbits, UGI nbitsperstage, UGI nbuckets
  ) : npushstages(npushstages), ntotalbits(ntotalbits), ntailbits(ntailbits), 
  nbitsperstage(nbitsperstage), nbuckets(nbuckets), 
  push_parts(npushstages+1), push_subparts(npushstages+1), 
  push_maps(npushstages+1), push_mapoffs(npushstages+1) {
    for(I i = 0; i < npushstages+1; ++i) {
      push_mapoffs[i].push_back(0);
    }
    pull_mapoffs.push_back(0);
  }
  
  void reset(UGI npushstages, UGI ntotalbits, UGI ntailbits, UGI nbitsperstage, UGI nbuckets) {
    this->npushstages = npushstages;
    this->ntotalbits = ntotalbits;
    this->ntailbits = ntailbits;
    this->nbitsperstage = nbitsperstage;
    this->nbuckets = nbuckets;
    
    input_parts.clear();
    push_parts.clear();
    push_subparts.clear();
    push_maps.clear();
    push_mapoffs.clear();
    pull_maps.clear();
    pull_mapoffs.clear();
    pull_parts.clear();
    
    push_parts.resize(npushstages+1);
    push_subparts.resize(npushstages+1);
    push_maps.resize(npushstages+1);
    push_mapoffs.resize(npushstages+1);
    
    for(UGI i = 0; i < npushstages+1; ++i) {
      push_mapoffs[i].push_back(0);
    }
    pull_mapoffs.push_back(0);
  }
  
  UGI get_max_push_buckets() const {
    return nbuckets;
  }
  
  void set_input_spaces(GI const count, ispace<GI> const * iparts) {
    for(GI i = 0; i < count; ++i)
      input_parts.push_back(iparts[i]);
  }
  
  GI get_num_input_spaces() const {
    return input_parts.size();
  }
  
  ispace<GI> const * get_input_spaces() const {
    return &input_parts[0];
  }
  
  UGI get_num_push_stages() const {
    return npushstages;
  }
  
  ispace<UGI> get_push_stage_space(UGI stage_idx) const {
    ispace<UGI> bounds(0,0);
    I nparts = push_parts[stage_idx].size();
    if(nparts > 0)
      bounds.set(push_parts[stage_idx][0].start, push_parts[stage_idx][nparts-1].end);
    return bounds;
  }
  
  UGI get_num_push_parts(UGI stage_idx) const {
    return push_parts[stage_idx].size();
  }
  
  ispace<UGI> const & get_push_part_space(UGI stage_idx, UGI part_idx) const {
    return push_parts[stage_idx][part_idx];
  }
  
  UGI append_push_part(UGI stage_idx, UGI part_size) {
    UGI numparts = push_parts[stage_idx].size();
    UGI start_idx = 0;
    UGI begin_subpart_idx = 0;
    
    if(numparts > 0) {
      start_idx = push_parts[stage_idx][numparts-1].end;
      begin_subpart_idx = push_subparts[stage_idx][numparts-1].end;
    }
    
    push_parts[stage_idx].push_back(
      ispace<UGI>(start_idx, start_idx + part_size)
    );
    
    push_subparts[stage_idx].push_back(
      ispace<UGI>(begin_subpart_idx, begin_subpart_idx)
    );
    
    push_mapoffs[stage_idx].push_back(push_mapoffs[stage_idx][numparts]);
    
    return push_parts[stage_idx].size()-1;
  }
  
  void set_subpart_range(
    UGI stage_idx, UGI part_idx, UGI begin_subpart_idx, UGI end_subpart_idx
  ) {
    push_subparts[stage_idx][part_idx].set(begin_subpart_idx, end_subpart_idx);
  }
  
  ispace<UGI> const & get_subpart_range(UGI stage_idx, UGI part_idx) const {
    return push_subparts[stage_idx][part_idx];
  }
  
  void set_push_map(UGI stage_idx, UGI part_idx, std::vector<DI> const & push_map) {
    UGI nparts = push_mapoffs[stage_idx].size()-1;
    UGI current_size = push_mapoffs[stage_idx][part_idx+1]-push_mapoffs[stage_idx][part_idx];
    UGI new_size = push_map.size();
    UGI startpos = push_mapoffs[stage_idx][part_idx];
    UGI endpos = push_mapoffs[stage_idx][part_idx+1];
    auto startiter = push_maps[stage_idx].begin()+startpos;
    auto enditer = push_maps[stage_idx].begin()+endpos;
    if(endpos != startpos) {
      push_maps[stage_idx].erase(startiter, enditer);
    }
    push_maps[stage_idx].insert(startiter, push_map.begin(), push_map.end());
    for(UGI i = part_idx; i < nparts; ++i) {
      push_mapoffs[stage_idx][i+1] += new_size;
      push_mapoffs[stage_idx][i+1] -= current_size;
    }
  }
  
  DI const * get_push_map(UGI stage_idx, UGI part_idx) const {
    return &push_maps[stage_idx][push_mapoffs[stage_idx][part_idx]];
  }
  
  size_t get_push_map_size(UGI stage_idx, UGI part_idx) const {
    return push_mapoffs[stage_idx][part_idx+1]-push_mapoffs[stage_idx][part_idx];
  }
  
  UGI get_num_pull_stages() const {
    return 1;
  }
  
  ispace<UGI> get_pull_stage_space() const {
    UGI npullparts = pull_parts.size();
    UGI start = 0, end = 0;
    if(npullparts > 0) {
      start = pull_parts[0].start;
      end = pull_parts[npullparts-1].end;
    }
    return ispace<UGI>(start, end);
  }
  
  UGI get_num_pull_parts() const {
    return pull_parts.size();
  }
  
  ispace<UGI> const & get_pull_part_space(UGI part_idx) const {
    return pull_parts[part_idx];
  }
  
  UGI append_pull_part(UGI part_size) {
    UGI npullparts = pull_parts.size();
    UGI start_idx = 0;
    if(npullparts > 0) {
      start_idx = pull_parts[npullparts-1].end;
    }
    pull_parts.push_back(ispace<UGI>(start_idx, start_idx+part_size));
    pull_mapoffs.push_back(pull_mapoffs[npullparts]);
    return pull_parts.size()-1;
  }
  
  void set_pull_map(UGI part_idx, std::vector<UGI> const & pull_map) {
    UGI nparts = pull_parts.size();
    UGI current_size = pull_mapoffs[part_idx+1]-pull_mapoffs[part_idx];
    UGI new_size = pull_map.size();
    UGI startpos = pull_mapoffs[part_idx];
    UGI endpos = pull_mapoffs[part_idx+1];
    auto startiter = pull_maps.begin()+startpos;
    auto enditer = pull_maps.begin()+endpos;
    if(endpos != startpos)
      pull_maps.erase(startiter, enditer);
    pull_maps.insert(startiter, pull_map.begin(), pull_map.end());
    for(UGI i = part_idx; i < nparts; ++i) {
      pull_mapoffs[i+1] += new_size;
      pull_mapoffs[i+1] -= current_size;
    }
  }
  
  UGI const * get_pull_map(UGI part_idx) const {
    return &pull_maps[pull_mapoffs[part_idx]];
  }
  
  size_t get_pull_map_size(UGI part_idx) const {
    return pull_mapoffs[part_idx+1]-pull_mapoffs[part_idx];
  }
  
  UGI get_max_buffer_size() const {
    UGI size = 0;
    for(UGI i = 0; i < npushstages+1; ++i) {
      UGI nparts = push_parts[i].size();
      size = std::max(size, push_parts[i][nparts-1].end-push_parts[i][0].start);
    }
    UGI npullparts = pull_parts.size();
    UGI pullsize = pull_parts[npullparts-1].end-pull_parts[0].start;
    return std::max(size, pullsize);
  }
};

template<typename GI, typename DI>
void record_distribution_sort(
  distribution_schedule<GI, DI> & schedule,
  typename std::make_unsigned<GI>::type const ntailbits,
  typename std::make_unsigned<GI>::type const nbitsperstage,
  typename std::make_unsigned<GI>::type const nstages,
  typename std::make_unsigned<GI>::type const stage_idx,
  typename std::make_unsigned<GI>::type const part_idx,
  typename std::make_unsigned<GI>::type const source_count,
  std::tuple<
    typename std::make_unsigned<GI>::type,
    typename std::make_unsigned<GI>::type,
    typename std::make_unsigned<GI>::type
  > * start,
  std::tuple<
    typename std::make_unsigned<GI>::type,
    typename std::make_unsigned<GI>::type,
    typename std::make_unsigned<GI>::type
  > * end
) {
  using UGI = typename std::make_unsigned<GI>::type;
  
  enum { BKT = 0, SRC, DST };
  
  if(stage_idx >= nstages) {
    // record permutation for each part in the output of final stage
    
    // renumber src
    UGI src_count = 0;
    UGI psrc = std::get<SRC>(*start);
    for(auto ptr = start; ptr < end; ++ptr) {
      UGI & src = std::get<SRC>(*ptr);
      if(psrc != src) {
        src_count++;
      }
      psrc = src = src_count;
    }
    
    // sort
    std::sort(start, end, [](auto const & lhs, auto const & rhs) {
        return std::get<DST>(lhs) < std::get<DST>(rhs);
      }
    );
    
    // check if dst is not repeated
    for(auto ptr = start+1; ptr < end; ++ptr) {
      UGI pdst = std::get<DST>(*(ptr-1));
      UGI dst = std::get<DST>(*ptr);
      if(dst == pdst) {
        printf("error: dst repeated\n");
      }
    }
    
    // form pull map
    std::vector<UGI> pullmap;
    for(auto ptr = start; ptr < end; ++ptr) {
      pullmap.push_back(std::get<SRC>(*ptr));
    }
    
    UGI pull_part_id = schedule.append_pull_part(pullmap.size());
    schedule.set_pull_map(pull_part_id, pullmap);
    
    return;
  }
  
  UGI const nbuckets = (UGI)1 << nbitsperstage;
  UGI const shift = ntailbits+nbitsperstage*(nstages-stage_idx-1);
  UGI const mask = (nbuckets-1) << shift;
  
  // calculate new bucket of each record
  // create distribution map
  std::vector<DI> distmap;
  {
    std::set<UGI> bucketset;
    
    UGI psrc = std::get<SRC>(*start);
    
    for(UGI i = 0; i < psrc; ++i)
      distmap.push_back(0);
    
    for(auto ptr = start; ptr < end; ++ptr) {
      UGI & bkt = std::get<BKT>(*ptr);
      UGI const & src = std::get<SRC>(*ptr);
      UGI const & dst = std::get<DST>(*ptr);
      bkt = (dst & mask) >> shift;
      if(psrc != src) {
        distmap.push_back(bucketset.size());
        for(auto b : bucketset) distmap.push_back(b);
        bucketset.clear();
        while(++psrc < src)
          distmap.push_back(0);
      }
      bucketset.insert(bkt);
    }
    distmap.push_back(bucketset.size());
    for(auto b : bucketset) distmap.push_back(b);
    
    while(++psrc < source_count)
      distmap.push_back(0);
  }
  
  std::sort(start, end);
  
  std::vector<UGI> bucket_boundaries(nbuckets+1, 0);
  std::vector<UGI> output_boundaries(nbuckets+1, 0);
  
  // determine bucket boundaries
  // determine output boundaries
  // renumber src
  {
    UGI psrc = std::get<SRC>(*start);
    UGI src_count = 0;
    
    UGI pbkt = std::get<BKT>(*start);
    auto bkt_start = start;
    auto ptr = start;
    
    while(ptr < end) {
      UGI const & bkt = std::get<BKT>(*ptr);
      UGI & src = std::get<SRC>(*ptr);
      
      if(pbkt != bkt) {
        bucket_boundaries[pbkt+1] = ptr-bkt_start;
        output_boundaries[pbkt+1] = src_count+1;
        
        bkt_start = ptr;
        pbkt = bkt;
        
        src_count = 0;
        psrc = src;
        src = src_count;
      } else {
        if(psrc != src) {
          src_count++;
          psrc = src;
          src = src_count;
        } else {
          src = src_count;
        }
      }
      
      ++ptr;
    }
    bucket_boundaries[pbkt+1] = ptr-bkt_start;
    output_boundaries[pbkt+1] = src_count+1;
    
    for(UGI b = 0; b < nbuckets; ++b) {
      bucket_boundaries[b+1] += bucket_boundaries[b];
      output_boundaries[b+1] += output_boundaries[b];
    }
  }
  
  std::vector<std::tuple<UGI, UGI, UGI>> subparts(nbuckets);
  UGI subpart_count = 0;
  for(UGI b = 0; b < nbuckets; ++b) {
    UGI subpart_start = output_boundaries[b];
    UGI subpart_end = output_boundaries[b+1];
    UGI subpart_size = subpart_end-subpart_start;
    if(subpart_size > 0) {
      UGI subpart_idx = schedule.append_push_part(stage_idx+1, subpart_size);
      subparts[subpart_count] = std::make_tuple(b, subpart_count, subpart_idx);
      subpart_count++;
    }
  }
  
  schedule.set_subpart_range(
    stage_idx, part_idx,
    std::get<2>(subparts[0]), std::get<2>(subparts[subpart_count-1])+1
  );
  
  schedule.set_push_map(stage_idx, part_idx, distmap);
  
  // recursively call on the subparts
  for(UGI i = 0; i < subpart_count; ++i) {
    UGI subpart_idx = std::get<2>(subparts[i]);
    UGI b = std::get<0>(subparts[i]);
    UGI bstart = bucket_boundaries[b];
    UGI bend = bucket_boundaries[b+1];
    
    record_distribution_sort(
      schedule,
      ntailbits,
      nbitsperstage,
      nstages,
      stage_idx+1,
      subpart_idx,
      output_boundaries[b+1]-output_boundaries[b],
      &start[bstart],
      &start[bend]
    );
  }
}

template<typename GI, typename DI>
void record_distribution_schedule(
  distribution_schedule<GI, DI> & schedule,
  typename std::make_unsigned<GI>::type const ntailbits,
  typename std::make_unsigned<GI>::type const nbitsperstage,
  GI const nparts,
  ispace<GI> const * parts,
  GI const outsize,
  GI const * out
) {
  using UGI = typename std::make_unsigned<GI>::type;
  
  UGI insize = 0;
  for(GI p = 0; p < nparts; ++p) {
    insize += parts[p].size;
  }
  
  // normalize out
  std::vector<UGI> normout(outsize);
  for(GI i = 0; i < outsize; ++i) {
    UGI idx = 0;
    for(GI p = 0; p < nparts; ++p) {
      if(parts[p].start <= out[i] && out[i] < parts[p].end) {
        idx += out[i]-parts[p].start;
        break;
      }
      idx += parts[p].size;
    }
    normout[i] = idx;
  }
  
  enum { BKT=0, SRC, DST };
  
  using info_type = std::tuple<UGI, UGI, UGI>;
  
  std::vector<info_type> info(outsize);
  for(GI i = 0; i < outsize; ++i) {
    std::get<BKT>(info[i]) = 0;
    std::get<SRC>(info[i]) = normout[i];
    std::get<DST>(info[i]) = i;
  }
  
  std::sort(info.begin(), info.end());
  
  UGI nbuckets = (UGI)1 << nbitsperstage;
  UGI ntotalbits = qsm::integer_log2((UGI)outsize)+1;
  UGI nstages = (
    ntotalbits > ntailbits ?
    (ntotalbits-ntailbits)/nbitsperstage + ((ntotalbits-ntailbits)%nbitsperstage ? 1 : 0)
    : 1
  );
  
  schedule.reset(nstages, ntotalbits, ntailbits, nbitsperstage, nbuckets);
  UGI part_idx = schedule.append_push_part(0, insize);
  schedule.set_input_spaces(nparts, parts);
  
  record_distribution_sort(
    schedule,
    ntailbits,
    nbitsperstage,
    nstages,
    0u,
    part_idx,
    insize,
    &info[0],
    &info[outsize]
  );
  
//#define REPLAY_DISTRIBUTION_SCHEDULE
#ifdef REPLAY_DISTRIBUTION_SCHEDULE
  ///////////////////////////////
  // replay the schedule to check
  ///////////////////////////////
  {
    // prepare input and output arrays
    auto maxsize = schedule.get_max_buffer_size();
    std::vector<GI> input(maxsize), output(maxsize);
    for(GI p = 0, idx = 0; p < nparts; ++p) {
      for(GI i = parts[p].start; i < parts[p].end; ++i) {
        input[idx++] = i;
      }
    }
    
    GI * iptr = &input[0];
    GI * optr = &output[0];
    GI * outptr[nbuckets];
    
    // distribution stages
    auto nstages = schedule.get_num_push_stages();
    for(decltype(nstages) s = 0; s < nstages; ++s) {
      // execute distribution schedule for each part in this stage
      auto nparts = schedule.get_num_push_parts(s);
      for(decltype(nparts) p = 0; p < nparts; ++p) {
        auto const & bounds = schedule.get_push_part_space(s, p);
        auto const & subparts = schedule.get_subpart_range(s, p);
        
        // prepare inputs and outputs
        auto inptr = &iptr[bounds.start];
        for(decltype(subparts.start) sp = subparts.start, idx = 0; sp < subparts.end; ++sp, ++idx) {
          auto spbounds = schedule.get_push_part_space(s+1, sp);
          outptr[idx] = &optr[spbounds.start];
        }
        
        // process input using distribution map to produce output
        auto distmap = schedule.get_push_map(s, p);
        for(decltype(bounds.size) i = 0; i < bounds.size; ++i) {
          auto nbkts = *distmap++;
          for(decltype(nbkts) j = 0; j < nbkts; ++j) {
            *(outptr[*distmap]) = inptr[i];
            outptr[*distmap]++;
            distmap++;
          }
        }
      }
      
      std::swap(iptr, optr);
    }
    
    // mapped gather
    auto nlastparts = schedule.get_num_push_parts(nstages);
    for(auto p = 0; p < nlastparts; ++p) {
      auto const & inpart = schedule.get_push_part_space(nstages, p);
      auto const & outpart = schedule.get_pull_part_space(p);
      auto const & gathermap = schedule.get_pull_map(p);
      
      GI * inptr = &iptr[inpart.start];
      GI * outptr = &optr[outpart.start];
      for(auto i = 0; i < outpart.size; ++i) {
        outptr[i] = inptr[gathermap[i]];
      }
    }
    
    // check errors
    GI err = 0;
    for(GI i = 0; i < outsize; ++i)
      if(out[i] != optr[i])
        err++;
    printf("Number of errors in distribution schedule = %u\n", (UGI)err);
  }
#endif // #ifdef REPLAY_DISTRIBUTION_SCHEDULE
}

template<typename GI, typename DI>
void record_distribution_schedule_for_seggregated_values(
  distribution_schedule<GI, DI> & schedule,
  typename std::make_unsigned<GI>::type const ntailbits,
  typename std::make_unsigned<GI>::type const nbitsperstage,
  GI const outsize,
  GI const * out
) {
  using UGI = typename std::make_unsigned<GI>::type;
  
  // normalize output
  std::set<GI> unique;
  for(GI i = 0; i < outsize; ++i) {
    unique.insert(out[i]);
  }
  
  UGI insize = unique.size();
  
  GI const smallest = *(unique.begin());
  GI const largest = *(unique.rbegin());
  
  std::vector<GI> o2nidx(largest-smallest+1, 0);
  {
    GI idx = 0;
    for(auto a : unique) {
      o2nidx[a-smallest] = idx++;
    }
  }
  
  std::vector<UGI> normout(outsize);
  for(GI i = 0; i < outsize; ++i) {
    normout[i] = o2nidx[out[i]-smallest];
  }
  
  // process normalized out stream
  enum { BKT=0, SRC, DST };
  
  using info_type = std::tuple<UGI, UGI, UGI>;
  
  std::vector<info_type> info(outsize);
  for(GI i = 0; i < outsize; ++i) {
    std::get<BKT>(info[i]) = 0;
    std::get<SRC>(info[i]) = normout[i];
    std::get<DST>(info[i]) = i;
  }
  
  std::sort(info.begin(), info.end());
  
  UGI nbuckets = (UGI)1 << nbitsperstage;
  UGI ntotalbits = qsm::integer_log2((UGI)outsize)+1;
  UGI nstages = (
    ntotalbits > ntailbits ?
    (ntotalbits-ntailbits)/nbitsperstage + ((ntotalbits-ntailbits)%nbitsperstage ? 1 : 0)
    : 1
  );
  //printf("nbuckets = %u, ntotalbits = %u, nstages = %u\n", nbuckets, ntotalbits, nstages);
  
  schedule.reset(nstages, ntotalbits, ntailbits, nbitsperstage, nbuckets);
  UGI part_idx = schedule.append_push_part(0, insize);
  GI nparts = 1;
  std::vector<ispace<GI>> parts = { ispace<GI>(0, insize) };
  
  schedule.set_input_spaces(nparts, &parts[0]);
  
  record_distribution_sort(
    schedule,
    ntailbits,
    nbitsperstage,
    nstages,
    0u,
    part_idx,
    insize,
    &info[0],
    &info[outsize]
  );
  
//#define REPLAY_DISTRIBUTION_SCHEDULE_FOR_SEGGREGATED_VALUES
#ifdef REPLAY_DISTRIBUTION_SCHEDULE_FOR_SEGGREGATED_VALUES
  ///////////////////////////////
  // replay the schedule to check
  ///////////////////////////////
  //{
  //  auto maxsize = schedule.get_max_buffer_size();
  //  std::vector<GI> input(maxsize), output(maxsize), tmp1(maxsize), tmp2(maxsize);
  //  GI idx = 0;
  //  for(auto a : unique)
  //    input[idx++] = a;
  //  
  //  execute_distribution_schedule(
  //    &output[0],
  //    &input[0],
  //    schedule,
  //    &tmp1[0],
  //    &tmp2[0]
  //  );
  //  
  //  GI err = 0;
  //  for(GI i = 0; i < outsize; ++i)
  //    if(out[i] != output[i])
  //      err++;
  //  printf("Number of errors in distribution schedule for seggregated values = %u\n", (UGI)err);
  //}
  {
    // prepare input and output arrays
    auto maxsize = schedule.get_max_buffer_size();
    std::vector<GI> input(maxsize), output(maxsize);
    {
      GI idx = 0;
      for(auto a : unique) {
        input[idx++] = a;
      }
    }
    
    GI * iptr = &input[0];
    GI * optr = &output[0];
    GI * outptr[nbuckets];
    
    // distribution stages
    auto nstages = schedule.get_num_push_stages();
    for(decltype(nstages) s = 0; s < nstages; ++s) {
      // execute distribution schedule for each part in this stage
      auto nparts = schedule.get_num_push_parts(s);
      for(decltype(nparts) p = 0; p < nparts; ++p) {
        auto const & bounds = schedule.get_push_part_space(s, p);
        auto const & subparts = schedule.get_subpart_range(s, p);
        
        // prepare inputs and outputs
        auto inptr = &iptr[bounds.start];
        for(decltype(subparts.start) sp = subparts.start, idx = 0; sp < subparts.end; ++sp, ++idx) {
          auto spbounds = schedule.get_push_part_space(s+1, sp);
          outptr[idx] = &optr[spbounds.start];
        }
        
        // process input using distribution map to produce output
        auto distmap = schedule.get_push_map(s, p);
        for(decltype(bounds.size) i = 0; i < bounds.size; ++i) {
          auto nbkts = *distmap++;
          for(decltype(nbkts) j = 0; j < nbkts; ++j) {
            *(outptr[*distmap]) = inptr[i];
            outptr[*distmap]++;
            distmap++;
          }
        }
      }
      
      std::swap(iptr, optr);
    }
    
    // mapped gather
    auto nlastparts = schedule.get_num_push_parts(nstages);
    for(auto p = 0; p < nlastparts; ++p) {
      auto const & inpart = schedule.get_push_part_space(nstages, p);
      auto const & outpart = schedule.get_pull_part_space(p);
      auto const & gathermap = schedule.get_pull_map(p);
      
      GI * inptr = &iptr[inpart.start];
      GI * outptr = &optr[outpart.start];
      for(auto i = 0; i < outpart.size; ++i) {
        outptr[i] = inptr[gathermap[i]];
      }
    }
    
    // check errors
    GI err = 0;
    for(GI i = 0; i < outsize; ++i)
      if(out[i] != optr[i])
        err++;
    printf("Number of errors in distribution schedule for seggregated values = %u\n", (UGI)err);
  }
#endif // #ifdef REPLAY_DISTRIBUTION_SCHEDULE_FOR_SEGGREGATED_VALUES
}

template<typename... T>
class distribution_buffer {
public:
  using tmp_type = std::tuple<T...>;
  using tmpptr_type = tmp_type *;
  using outptr_type = std::tuple<T *...>;
  using const_outptr_type = std::tuple<T const *...>;
  using inptr_type = std::tuple<T const *...>;
  
private:
  tmpptr_type * m_tmp_buffer[2];
  outptr_type * m_out_buffer;
  
  distribution_buffer(distribution_buffer<T...> const &) = delete;
  distribution_buffer<T...> & operator=(distribution_buffer<T...> const &) = delete;
  
public:
  distribution_buffer() : m_tmp_buffer{nullptr, nullptr}, m_out_buffer(nullptr) {
  }
  
  ~distribution_buffer() {
    deallocate();
  }
  
  tmpptr_type tmp(int i, int p) const {
    return m_tmp_buffer[i][p];
  }
  
  outptr_type const & operator[](int p) {
    return m_out_buffer[p];
  }
  
  const_outptr_type operator[](int p) const {
    return m_out_buffer[p];
  }
  
  template<typename GI, typename DI>
  void allocate(int nparts, distribution_schedule<GI, DI> const * schedules) {
    deallocate();
    
    std::vector<GI> bufferoffs(nparts+1), outoffs(nparts+1);
    bufferoffs[0] = 0;
    outoffs[0] = 0;
    for(GI p = 0; p < nparts; ++p) {
      bufferoffs[p+1] = bufferoffs[p] + schedules[p].get_max_buffer_size();
      outoffs[p+1] = outoffs[p] + schedules[p].get_pull_stage_space().size;
    }
    
    m_tmp_buffer[0] = new std::tuple<T...>*[nparts];
    m_tmp_buffer[1] = new std::tuple<T...>*[nparts];
    m_out_buffer = new std::tuple<T * ...>[nparts];
    std::tuple<T...> * buffer1_mem = new std::tuple<T...>[bufferoffs[nparts]];
    std::tuple<T...> * buffer2_mem = new std::tuple<T...>[bufferoffs[nparts]];
    std::tuple<T * ...> out_mem = { new T[outoffs[nparts]] ... };
    for(GI p = 0; p < nparts; ++p) {
      m_tmp_buffer[0][p] = &buffer1_mem[bufferoffs[p]];
      m_tmp_buffer[1][p] = &buffer2_mem[bufferoffs[p]];
      GI outoff = outoffs[p];
      m_out_buffer[p] = qsm::apply_with_tuple_result(
        [outoff](auto & elm) {
          return &elm[outoff];
        },
        out_mem
      );
    }
  }
  
  void deallocate() {
    if(m_tmp_buffer[0]) {
      delete[] m_tmp_buffer[0][0];
      delete[] m_tmp_buffer[0];
    }
    
    if(m_tmp_buffer[1]) {
      delete[] m_tmp_buffer[1][0];
      delete[] m_tmp_buffer[1];
    }
    
    if(m_out_buffer) {
      qsm::apply_with_tuple_result(
        [](auto & elm) { delete[] elm; },
        m_out_buffer[0]
      );
      delete[] m_out_buffer;
    }
    
    m_tmp_buffer[0] = nullptr;
    m_tmp_buffer[1] = nullptr;
    m_out_buffer = nullptr;
  }
};

template<typename T>
class distribution_buffer<T> {
public:
  using tmp_type = T;
  using tmpptr_type = tmp_type *;
  using outptr_type = T *;
  using inptr_type = T const *;
  
private:
  tmpptr_type * m_tmp_buffer[2];
  outptr_type * m_out_buffer;
  
public:
  distribution_buffer() : m_tmp_buffer{nullptr, nullptr}, m_out_buffer(nullptr) {
  }
  
  tmpptr_type tmp(int i, int p) const {
    return m_tmp_buffer[i][p];
  }
  
  outptr_type operator[](int p) const {
    return m_out_buffer[p];
  }
  
  template<typename GI, typename DI>
  void allocate(int nparts, distribution_schedule<GI, DI> const * schedules) {
    std::vector<GI> bufferoffs(nparts+1), outoffs(nparts+1);
    bufferoffs[0] = 0;
    outoffs[0] = 0;
    for(GI p = 0; p < nparts; ++p) {
      bufferoffs[p+1] = bufferoffs[p] + schedules[p].get_max_buffer_size();
      outoffs[p+1] = outoffs[p] + schedules[p].get_pull_stage_space().size;
    }
    
    m_tmp_buffer[0] = new T*[nparts];
    m_tmp_buffer[1] = new T*[nparts];
    m_out_buffer = new T*[nparts];
    T * buffer1_mem = new T[bufferoffs[nparts]];
    T * buffer2_mem = new T[bufferoffs[nparts]];
    T * out_mem = new T[outoffs[nparts]];
    for(GI p = 0; p < nparts; ++p) {
      m_tmp_buffer[0][p] = &buffer1_mem[bufferoffs[p]];
      m_tmp_buffer[1][p] = &buffer2_mem[bufferoffs[p]];
      m_out_buffer[p] = &out_mem[outoffs[p]];
    }
  }
  
  void deallocate() {
    if(m_tmp_buffer[0]) {
      delete[] m_tmp_buffer[0][0];
      delete[] m_tmp_buffer[0];
    }
    
    if(m_tmp_buffer[1]) {
      delete[] m_tmp_buffer[1][0];
      delete[] m_tmp_buffer[1];
    }
    
    if(m_out_buffer) {
      delete[] m_out_buffer[0];
      delete[] m_out_buffer;
    }
    
    m_tmp_buffer[0] = nullptr;
    m_tmp_buffer[1] = nullptr;
    m_out_buffer = nullptr;
  }
};

template<typename GI, typename... T, size_t... I>
inline std::tuple<T...> gather_details(
  qsm::seq_cgrqueue<GI, T...> const & q,
  std::index_sequence<I...>
) {
  return std::make_tuple(
    (q.template get<I>()) ...
  );
}

template<typename GI, typename... T>
inline std::tuple<T...> gather(qsm::seq_cgrqueue<GI, T...> const & q) {
  #pragma forceinline
  return gather_details(q, std::make_index_sequence<sizeof...(T)>{});
}

template<typename GI, typename... T, size_t... I>
inline void scatter_details(
  std::tuple<T *...> const & ptr,
  int const idx,
  std::tuple<T...> const & in,
  std::index_sequence<I...>
) {
  auto r = std::make_tuple( (std::get<I>(ptr)[idx] = std::get<I>(in)) ... );
}

template<typename GI, typename... T>
inline void scatter(
  std::tuple<T *...> const & ptr,
  int const idx,
  std::tuple<T...> const & in
) {
  #pragma forceinline
  scatter_details<GI, T...>(ptr, idx, in, std::make_index_sequence<sizeof...(T)>{});
}

template<typename T, typename GI, typename DI>
void execute_distribution_schedule(
  T * output,
  T const * input,
  distribution_schedule<GI, DI> const & schedule,
  T * buffer1,
  T * buffer2
) {
  using UGI = typename distribution_schedule<GI, DI>::uindex_type;
  
  GI const ninspc = schedule.get_num_input_spaces();
  ispace<GI> const * inspc = schedule.get_input_spaces();
  
  UGI const npushstages = schedule.get_num_push_stages();
  
  T const * iptr;
  T * optr;
  T * outptr[schedule.get_max_push_buckets()];
  
  // push phase
  if(npushstages > 0) {
    // execute first push stage
    {
      optr = buffer1;
      
      ispace<UGI> const & subparts = schedule.get_subpart_range(0, 0);
      for(UGI sp = subparts.start, idx = 0; sp < subparts.end; ++sp, ++idx) {
        ispace<UGI> const & subpart_space = schedule.get_push_part_space(1, sp);
        outptr[idx] = &optr[subpart_space.start];
      }
      
      DI const * push_map = schedule.get_push_map(0, 0);
      
      for(GI i = 0; i < ninspc; ++i) {
        T const * inptr = &input[inspc[i].start];
        
        for(GI j = 0; j < inspc[i].size; ++j) {
          DI nbkts = *push_map++;
          for(DI k = 0; k < nbkts; ++k) {
            *(outptr[*push_map]) = inptr[j];
            outptr[*push_map]++;
            push_map++;
          }
        }
      }
      
      iptr = optr;
      optr = buffer2;
    }
    
    // execute remaining push stages
    for(UGI s = 1; s < npushstages; ++s) {
      UGI const nparts = schedule.get_num_push_parts(s);
      for(UGI p = 0; p < nparts; ++p) {
        ispace<UGI> const & part_space = schedule.get_push_part_space(s, p);
        ispace<UGI> const & subpart_range = schedule.get_subpart_range(s, p);
        
        T const * inptr = &iptr[part_space.start];
        for(UGI sp = subpart_range.start, idx = 0; sp < subpart_range.end; ++sp, ++idx) {
          ispace<UGI> const & subpart_space = schedule.get_push_part_space(s+1, sp);
          outptr[idx] = &optr[subpart_space.start];
        }
        
        DI const * push_map = schedule.get_push_map(s, p);
        
        for(UGI i = 0; i < part_space.size; ++i) {
          DI nbkts = *push_map++;
          for(DI j = 0; j < nbkts; ++j) {
            *(outptr[*push_map]) = inptr[i];
            outptr[*push_map]++;
            push_map++;
          }
        }
      }
      
      T * tptr = const_cast<T *>(iptr);
      iptr = optr;
      optr = tptr;
    }
  } else {
    optr = buffer1;
    
    for(GI i = 0, idx = 0; i < ninspc; ++i) {
      T const * inptr = &input[inspc[i].start];
      
      T * outptr = &optr[idx];
      for(GI j = 0; j < inspc[i].size; ++j)
        optr[j] = inptr[j];
      idx += inspc[i].size;
    }
    
    iptr = optr;
    optr = buffer2;
  }
  
  // pull phase
  UGI const npullparts = schedule.get_num_pull_parts();
  for(UGI p = 0; p < npullparts; ++p) {
    ispace<UGI> const & inspace = schedule.get_push_part_space(npushstages, p);
    ispace<UGI> const & outspace = schedule.get_pull_part_space(p);
    
    UGI const * pull_map = schedule.get_pull_map(p);
    T const * inptr = &iptr[inspace.start];
    T * outptr = &output[outspace.start];
    
    for(UGI i = 0; i < outspace.size; ++i) {
      outptr[i] = inptr[pull_map[i]];
    }
  }
}

template<typename GI, typename DI, typename... T>
void execute_distribution_schedule_in_pack(
  std::tuple<T * ...> const & output,
  std::tuple<T const * ...> const & input,
  distribution_schedule<GI, DI> const & schedule,
  std::tuple<T ...> * buffer1,
  std::tuple<T ...> * buffer2
) {
  using VT = std::tuple<T...>;
  using UGI = typename distribution_schedule<GI, DI>::uindex_type;
  
  GI const ninspc = schedule.get_num_input_spaces();
  ispace<GI> const * inspc = schedule.get_input_spaces();
  
  UGI const npushstages = schedule.get_num_push_stages();
  
  VT * iptr = nullptr;
  VT * optr = nullptr;
  
  // push phase
  if(npushstages > 0) {
    VT * outptr[schedule.get_max_push_buckets()];
    
    // execute first push stage
    {
      optr = buffer1;
      
      ispace<UGI> const & subparts = schedule.get_subpart_range(0, 0);
      for(UGI sp = subparts.start, idx = 0; sp < subparts.end; ++sp, ++idx) {
        ispace<UGI> const & subpart_space = schedule.get_push_part_space(1, sp);
        outptr[idx] = &optr[subpart_space.start];
      }
      
      DI const * push_map = schedule.get_push_map(0, 0);
      
      for(GI i = 0; i < ninspc; ++i) {
        std::tuple<T const * ...> inptr = qsm::pointers_at(input, inspc[i].start);
        
        for(GI j = 0; j < inspc[i].size; ++j) {
          DI nbkts = *push_map++;
          for(DI k = 0; k < nbkts; ++k) {
            #pragma forceinline
            qsm::assign(*(outptr[*push_map]), inptr, j);
            outptr[*push_map]++;
            push_map++;
          }
        }
      }
      
      iptr = optr;
      optr = buffer2;
    }
    
    // execute remaining push stages
    for(UGI s = 1; s < npushstages; ++s) {
      UGI const nparts = schedule.get_num_push_parts(s);
      for(UGI p = 0; p < nparts; ++p) {
        ispace<UGI> const & part_space = schedule.get_push_part_space(s, p);
        ispace<UGI> const & subpart_range = schedule.get_subpart_range(s, p);
        
        VT const * inptr = &iptr[part_space.start];
        for(UGI sp = subpart_range.start, idx = 0; sp < subpart_range.end; ++sp, ++idx) {
          ispace<UGI> const & subpart_space = schedule.get_push_part_space(s+1, sp);
          outptr[idx] = &optr[subpart_space.start];
        }
        
        DI const * push_map = schedule.get_push_map(s, p);
        
        for(UGI i = 0; i < part_space.size; ++i) {
          DI nbkts = *push_map++;
          for(DI j = 0; j < nbkts; ++j) {
            *(outptr[*push_map]) = inptr[i];
            outptr[*push_map]++;
            push_map++;
          }
        }
      }
      
      VT * tptr = const_cast<VT *>(iptr);
      iptr = optr;
      optr = tptr;
    }
  } else {
    optr = buffer1;
    
    for(GI i = 0, idx = 0; i < ninspc; ++i) {
      std::tuple<T const * ...> inptr = qsm::pointers_at(input, inspc[i].start);
      
      VT * outptr = &optr[idx];
      for(GI j = 0; j < inspc[i].size; ++j) {
        #pragma forceinline
        qsm::assign(optr[j], inptr, j);
      }
      idx += inspc[i].size;
    }
    
    iptr = optr;
    optr = buffer2;
  }
  
  // pull phase
  UGI const npullparts = schedule.get_num_pull_parts();
  for(UGI p = 0; p < npullparts; ++p) {
    ispace<UGI> const & inspace = schedule.get_push_part_space(npushstages, p);
    ispace<UGI> const & outspace = schedule.get_pull_part_space(p);
    
    UGI const * pull_map = schedule.get_pull_map(p);
    VT const * inptr = &iptr[inspace.start];
    std::tuple<T * ...> outptr = qsm::pointers_at(output, outspace.start);
    
    for(UGI i = 0; i < outspace.size; ++i) {
      #pragma forceinline
      qsm::assign(outptr, i, inptr[pull_map[i]]);
    }
  }
}

template<typename GI, typename DI, typename... T>
void execute_distribution_schedule_in_pack(
  std::tuple<T * ...> const & output,
  std::tuple<T ...> const * input,
  distribution_schedule<GI, DI> const & schedule,
  std::tuple<T ...> * buffer1,
  std::tuple<T ...> * buffer2
) {
  using VT = std::tuple<T...>;
  using UGI = typename distribution_schedule<GI, DI>::uindex_type;
  
  GI const ninspc = schedule.get_num_input_spaces();
  ispace<GI> const * inspc = schedule.get_input_spaces();
  
  UGI const npushstages = schedule.get_num_push_stages();
  
  VT * iptr = nullptr;
  VT * optr = nullptr;
  
  // push phase
  if(npushstages > 0) {
    VT * outptr[schedule.get_max_push_buckets()];
    
    // execute first push stage
    {
      optr = buffer1;
      
      ispace<UGI> const & subparts = schedule.get_subpart_range(0, 0);
      for(UGI sp = subparts.start, idx = 0; sp < subparts.end; ++sp, ++idx) {
        ispace<UGI> const & subpart_space = schedule.get_push_part_space(1, sp);
        outptr[idx] = &optr[subpart_space.start];
      }
      
      DI const * push_map = schedule.get_push_map(0, 0);
      
      for(GI i = 0; i < ninspc; ++i) {
        std::tuple<T ...> const * inptr = &input[inspc[i].start];
        
        for(GI j = 0; j < inspc[i].size; ++j) {
          DI nbkts = *push_map++;
          for(DI k = 0; k < nbkts; ++k) {
            //#pragma forceinline
            //qsm::assign(*(outptr[*push_map]), inptr, j);
            *(outptr[*push_map]) = inptr[j];
            outptr[*push_map]++;
            push_map++;
          }
        }
      }
      
      iptr = optr;
      optr = buffer2;
    }
    
    // execute remaining push stages
    for(UGI s = 1; s < npushstages; ++s) {
      UGI const nparts = schedule.get_num_push_parts(s);
      for(UGI p = 0; p < nparts; ++p) {
        ispace<UGI> const & part_space = schedule.get_push_part_space(s, p);
        ispace<UGI> const & subpart_range = schedule.get_subpart_range(s, p);
        
        VT const * inptr = &iptr[part_space.start];
        for(UGI sp = subpart_range.start, idx = 0; sp < subpart_range.end; ++sp, ++idx) {
          ispace<UGI> const & subpart_space = schedule.get_push_part_space(s+1, sp);
          outptr[idx] = &optr[subpart_space.start];
        }
        
        DI const * push_map = schedule.get_push_map(s, p);
        
        for(UGI i = 0; i < part_space.size; ++i) {
          DI nbkts = *push_map++;
          for(DI j = 0; j < nbkts; ++j) {
            *(outptr[*push_map]) = inptr[i];
            outptr[*push_map]++;
            push_map++;
          }
        }
      }
      
      VT * tptr = const_cast<VT *>(iptr);
      iptr = optr;
      optr = tptr;
    }
  } else {
    optr = buffer1;
    
    for(GI i = 0, idx = 0; i < ninspc; ++i) {
      std::tuple<T ...> const * inptr = &input[inspc[i].start];
      
      VT * outptr = &optr[idx];
      for(GI j = 0; j < inspc[i].size; ++j) {
        //#pragma forceinline
        //qsm::assign(optr[j], inptr, j);
        optr[j] = inptr[j];
      }
      idx += inspc[i].size;
    }
    
    iptr = optr;
    optr = buffer2;
  }
  
  // pull phase
  UGI const npullparts = schedule.get_num_pull_parts();
  for(UGI p = 0; p < npullparts; ++p) {
    ispace<UGI> const & inspace = schedule.get_push_part_space(npushstages, p);
    ispace<UGI> const & outspace = schedule.get_pull_part_space(p);
    
    UGI const * pull_map = schedule.get_pull_map(p);
    VT const * inptr = &iptr[inspace.start];
    std::tuple<T * ...> outptr = qsm::pointers_at(output, outspace.start);
    
    for(UGI i = 0; i < outspace.size; ++i) {
      #pragma forceinline
      qsm::assign(outptr, i, inptr[pull_map[i]]);
    }
  }
}

template<int BufferSize = 64, typename GI, typename DI, typename T>
void execute_distribution_schedule_using_qsm0(
  T * output,
  T const * input,
  distribution_schedule<GI, DI> const & schedule,
  T * buffer1,
  T * buffer2
) {
  using UGI = typename distribution_schedule<GI, DI>::uindex_type;
  
  GI const ninspc = schedule.get_num_input_spaces();
  ispace<GI> const * inspc = schedule.get_input_spaces();
  
  UGI const npushstages = schedule.get_num_push_stages();
  
  T const * iptr;
  T * optr;
  
  // push phase
  if(npushstages > 0) {
    qsm::seq_cgrqueue<GI, T> inq;
    qsm::cgwqueue<GI, T, BufferSize> outq[schedule.get_max_push_buckets()];
    T * bufptrs[schedule.get_max_push_buckets()];
    int bsz[schedule.get_max_push_buckets()];
    int bufidx[schedule.get_max_push_buckets()];
    
    // execute first push stage
    {
      optr = buffer1;
      
      ispace<UGI> const & subparts = schedule.get_subpart_range(0, 0);
      
      for(UGI sp = subparts.start, idx = 0; sp < subparts.end; ++sp, ++idx) {
        ispace<UGI> const & subpart_space = schedule.get_push_part_space(1, sp);
        outq[idx].configure(&optr[subpart_space.start]);
      }
      
      DI const * push_map = schedule.get_push_map(0, 0);
      
      for(GI i = 0; i < ninspc; ++i) {
        T const * inptr = &input[inspc[i].start];
        
        inq.configure(inptr);
        
        for(GI j = 0; j < inspc[i].size; j+=BufferSize) {
          for(int k = 0; k < subparts.size; ++k) {
            bufptrs[k] = outq[k].get(&bsz[k]);
            bufidx[k] = 0;
          }
          
          int niter = std::min(inspc[i].size-j, BufferSize);
          for(int l = 0; l < niter; ++l) {
            auto val = inq.pop();
            DI nbkts = *push_map++;
            for(DI k = 0; k < nbkts; ++k) {
              bufptrs[*push_map][bufidx[*push_map]] = val; //inptr[j+l];
              bufidx[*push_map]++;
              push_map++;
            }
          }
          
          for(int k = 0; k < subparts.size; ++k) {
            outq[k].flush(bufidx[k]);
          }
        }
      }
      
      for(UGI i = 0; i < subparts.size; ++i)
        outq[i].deactivate();
      
      iptr = optr;
      optr = buffer2;
      
      _mm_mfence();
    }
    
    // execute remaining push stages
    for(UGI s = 1; s < npushstages; ++s) {
      UGI const nparts = schedule.get_num_push_parts(s);
      for(UGI p = 0; p < nparts; ++p) {
        ispace<UGI> const & part_space = schedule.get_push_part_space(s, p);
        ispace<UGI> const & subparts = schedule.get_subpart_range(s, p);
        
        T const * inptr = &iptr[part_space.start];
        inq.configure(inptr);
        
        for(UGI sp = subparts.start, idx = 0; sp < subparts.end; ++sp, ++idx) {
          ispace<UGI> const & subpart_space = schedule.get_push_part_space(s+1, sp);
          outq[idx].configure(&optr[subpart_space.start]);
        }
        
        DI const * push_map = schedule.get_push_map(s, p);
        
        for(UGI i = 0; i < part_space.size; i+=BufferSize) {
          for(UGI k = 0; k < subparts.size; ++k) {
            bufptrs[k] = outq[k].get(&bsz[k]);
            bufidx[k] = 0;
          }
          
          int niter = std::min((int)(part_space.size-i), BufferSize);
          
          for(int j = 0; j < niter; ++j) {
            auto val = inq.pop();
            DI nbkts = *push_map++;
            for(DI k = 0; k < nbkts; ++k) {
              bufptrs[*push_map][bufidx[*push_map]] = val; //inptr[i+j];
              bufidx[*push_map]++;
              push_map++;
            }
          }
          
          for(UGI k = 0; k < subparts.size; ++k) {
            outq[k].flush(bufidx[k]);
          }
        }
        
        for(UGI k = 0; k < subparts.size; ++k) {
          outq[k].deactivate();
        }
      }
      
      T * tptr = const_cast<T *>(iptr);
      iptr = optr;
      optr = tptr;
      
      _mm_mfence();
    }
  } else {
    optr = buffer1;
    
    qsm::seq_cgrqueue<GI, T> inq;
    
    qsm::cgwqueue<GI, T, BufferSize> outq;
    T * bufptr;
    int bsz;
    
    outq.configure(optr);
    
    for(GI i = 0; i < ninspc; ++i) {
      T const * inptr = &input[inspc[i].start];
      
      inq.configure(inptr);
      
      for(GI j = 0; j < inspc[i].size; j+=BufferSize) {
        bufptr = outq.get(&bsz);
        int niter = std::min(inspc[i].size-j, BufferSize);
        for(int l = 0; l < niter; ++l) {
          bufptr[l] = inq.pop(); //inptr[j+l];
        }
        outq.flush(niter);
      }
      
      outq.deactivate();
    }
    
    iptr = optr;
    optr = buffer2;
    
    _mm_mfence();
  }
  
  // pull phase
  {
    UGI const npullparts = schedule.get_num_pull_parts();
    
    qsm::cgwqueue<GI, T, BufferSize> outq;
    T * bufptr;
    int bsz;
    
    outq.configure(output);
    
    for(UGI p = 0; p < npullparts; ++p) {
      ispace<UGI> const & inspace = schedule.get_push_part_space(npushstages, p);
      ispace<UGI> const & outspace = schedule.get_pull_part_space(p);
      
      UGI const * pull_map = schedule.get_pull_map(p);
      T const * inptr = &iptr[inspace.start];
      //T * outptr = &output[outspace.start];
      
      for(UGI i = 0; i < outspace.size; i+=BufferSize) {
        bufptr = outq.get();
        int niter = std::min((int)(outspace.size-i), BufferSize);
        #pragma loop_count max(BufferSize)
        for(int j = 0; j < niter; ++j) {
          bufptr[j] = inptr[pull_map[i+j]];
        }
        outq.flush(niter);
      }
    }
    
    outq.deactivate();
  }
}

template<int BufferSize = 64, typename GI, typename DI, typename... T>
void execute_distribution_schedule_in_pack_using_qsm0(
  std::tuple<T * ...> const & output,
  std::tuple<T const * ...> const & input,
  distribution_schedule<GI, DI> const & schedule,
  std::tuple<T...> * buffer1,
  std::tuple<T...> * buffer2
) {
  using VT = std::tuple<T...>;
  using UGI = typename distribution_schedule<GI, DI>::uindex_type;
  
  GI const ninspc = schedule.get_num_input_spaces();
  ispace<GI> const * inspc = schedule.get_input_spaces();
  
  UGI const npushstages = schedule.get_num_push_stages();
  
  VT * iptr = nullptr;
  VT * optr = nullptr;
  
  // push phase
  if(npushstages > 0) {
    qsm::cgwqueue<GI, VT, BufferSize> outq[schedule.get_max_push_buckets()];
    VT * bufptrs[schedule.get_max_push_buckets()];
    int bufidx[schedule.get_max_push_buckets()];
    
    // execute first push stage
    {
      optr = buffer1;
      
      qsm::seq_cgrqueue<GI, T...> inq;
      
      ispace<UGI> const & subparts = schedule.get_subpart_range(0, 0);
      
      for(UGI sp = subparts.start, idx = 0; sp < subparts.end; ++sp, ++idx) {
        ispace<UGI> const & subpart_space = schedule.get_push_part_space(1, sp);
        outq[idx].configure(&optr[subpart_space.start]);
      }
      
      DI const * push_map = schedule.get_push_map(0, 0);
      
      for(GI i = 0; i < ninspc; ++i) {
        inq.configure(inspc[i].start, input);
        
        for(GI j = 0; j < inspc[i].size; j+=BufferSize) {
          #pragma novector
          for(int k = 0; k < subparts.size; ++k) {
            bufptrs[k] = outq[k].get();
            bufidx[k] = 0;
          }
          
          int niter = std::min(inspc[i].size-j, BufferSize);
          
          #pragma loop_count max(BufferSize)
          for(int l = 0; l < niter; ++l) {
            DI nbkts = *push_map++;
            for(DI k = 0; k < nbkts; ++k) {
              #pragma forceinline
              bufptrs[*push_map][bufidx[*push_map]] = gather(inq);
              bufidx[*push_map]++;
              push_map++;
            }
            
            inq.advance();
          }
          
          for(int k = 0; k < subparts.size; ++k) {
            #pragma forceinline recursive
            outq[k].flush(bufidx[k]);
          }
        }
      }
      
      for(UGI i = 0; i < subparts.size; ++i)
        outq[i].deactivate();
      
      _mm_mfence();
    }
    
    iptr = buffer1;
    optr = buffer2;
    
    // execute remaining push stages
    qsm::seq_cgrqueue<GI, VT> inq;
    for(UGI s = 1; s < npushstages; ++s) {
      UGI const nparts = schedule.get_num_push_parts(s);
      for(UGI p = 0; p < nparts; ++p) {
        ispace<UGI> const & part_space = schedule.get_push_part_space(s, p);
        ispace<UGI> const & subparts = schedule.get_subpart_range(s, p);
        
        VT const * inptr = &iptr[part_space.start];
        inq.configure(inptr);
        
        for(UGI sp = subparts.start, idx = 0; sp < subparts.end; ++sp, ++idx) {
          ispace<UGI> const & subpart_space = schedule.get_push_part_space(s+1, sp);
          outq[idx].configure(&optr[subpart_space.start]);
        }
        
        DI const * push_map = schedule.get_push_map(s, p);
        
        for(UGI i = 0; i < part_space.size; i+=BufferSize) {
          #pragma novector
          for(UGI k = 0; k < subparts.size; ++k) {
            bufptrs[k] = outq[k].get();
            bufidx[k] = 0;
          }
          
          int niter = std::min((int)(part_space.size-i), BufferSize);
          
          #pragma loop_count max(BufferSize)
          for(int j = 0; j < niter; ++j) {
            auto val = inq.pop();
            DI nbkts = *push_map++;
            for(DI k = 0; k < nbkts; ++k) {
              bufptrs[*push_map][bufidx[*push_map]] = val; //inptr[i+j];
              bufidx[*push_map]++;
              push_map++;
            }
          }
          
          for(UGI k = 0; k < subparts.size; ++k) {
            #pragma forceinline recursive
            outq[k].flush(bufidx[k]);
          }
        }
        
        for(UGI k = 0; k < subparts.size; ++k) {
          outq[k].deactivate();
        }
      }
      
      VT * tptr = const_cast<VT *>(iptr);
      iptr = optr;
      optr = tptr;
      
      _mm_mfence();
    }
  } else {
    optr = buffer1;
    
    qsm::seq_cgrqueue<GI, T...> inq;
    
    qsm::cgwqueue<GI, VT, BufferSize> outq;
    VT * bufptr;
    
    outq.configure(optr);
    
    for(GI i = 0; i < ninspc; ++i) {
      inq.configure(inspc[i].start, input);
      
      for(GI j = 0; j < inspc[i].size; j+=BufferSize) {
        bufptr = outq.get();
        int niter = std::min(inspc[i].size-j, BufferSize);
        for(int l = 0; l < niter; ++l) {
          bufptr[l] = gather(inq);//inq.pop();
          inq.advance();
        }
        
        outq.flush(niter);
      }
      
      outq.deactivate();
    }
    
    iptr = optr;
    
    _mm_mfence();
  }
  
  // pull phase
  {
    UGI const npullparts = schedule.get_num_pull_parts();
    
    qsm::cgwqueue_pack<BufferSize, GI, T...> outq;
    std::tuple<T *...> bufptr;
    
    outq.configure(output);
    
    for(UGI p = 0; p < npullparts; ++p) {
      ispace<UGI> const & inspace = schedule.get_push_part_space(npushstages, p);
      ispace<UGI> const & outspace = schedule.get_pull_part_space(p);
      
      UGI const * pull_map = schedule.get_pull_map(p);
      VT const * inptr = &iptr[inspace.start];
      
      for(UGI i = 0; i < outspace.size; i+=BufferSize) {
        bufptr = outq.get();
        
        int niter = std::min((int)(outspace.size-i), BufferSize);
        
        #pragma loop_count max(BufferSize)
        for(int j = 0; j < niter; ++j) {
          #pragma forceinline
          scatter<GI, T...>(bufptr, j, inptr[pull_map[i+j]]);
        }
        
        outq.flush(niter);
      }
    }
    
    outq.deactivate();
  }
}

template<int BufferSize = 64, typename GI, typename DI, typename... T>
void execute_distribution_schedule_in_pack_using_qsm0(
  std::tuple<T * ...> const & output,
  std::tuple<T ...> const * input,
  distribution_schedule<GI, DI> const & schedule,
  std::tuple<T...> * buffer1,
  std::tuple<T...> * buffer2
) {
  using VT = std::tuple<T...>;
  using UGI = typename distribution_schedule<GI, DI>::uindex_type;
  
  GI const ninspc = schedule.get_num_input_spaces();
  ispace<GI> const * inspc = schedule.get_input_spaces();
  
  UGI const npushstages = schedule.get_num_push_stages();
  
  VT * iptr = nullptr;
  VT * optr = nullptr;
  
  // push phase
  if(npushstages > 0) {
    qsm::cgwqueue<GI, VT, BufferSize> outq[schedule.get_max_push_buckets()];
    VT * bufptrs[schedule.get_max_push_buckets()];
    int bufidx[schedule.get_max_push_buckets()];
    
    // execute first push stage
    {
      optr = buffer1;
      
      qsm::seq_cgrqueue<GI, VT> inq;
      
      ispace<UGI> const & subparts = schedule.get_subpart_range(0, 0);
      
      for(UGI sp = subparts.start, idx = 0; sp < subparts.end; ++sp, ++idx) {
        ispace<UGI> const & subpart_space = schedule.get_push_part_space(1, sp);
        outq[idx].configure(&optr[subpart_space.start]);
      }
      
      DI const * push_map = schedule.get_push_map(0, 0);
      
      for(GI i = 0; i < ninspc; ++i) {
        inq.configure(inspc[i].start, input);
        
        for(GI j = 0; j < inspc[i].size; j+=BufferSize) {
          #pragma novector
          for(int k = 0; k < subparts.size; ++k) {
            bufptrs[k] = outq[k].get();
            bufidx[k] = 0;
          }
          
          int niter = std::min(inspc[i].size-j, BufferSize);
          
          #pragma loop_count max(BufferSize)
          for(int l = 0; l < niter; ++l) {
            auto val = inq.pop();
            DI nbkts = *push_map++;
            for(DI k = 0; k < nbkts; ++k) {
              #pragma forceinline
              bufptrs[*push_map][bufidx[*push_map]] = val;//gather(inq);
              bufidx[*push_map]++;
              push_map++;
            }
            
            //inq.advance();
          }
          
          for(int k = 0; k < subparts.size; ++k) {
            #pragma forceinline recursive
            outq[k].flush(bufidx[k]);
          }
        }
      }
      
      for(UGI i = 0; i < subparts.size; ++i)
        outq[i].deactivate();
      
      _mm_mfence();
    }
    
    iptr = buffer1;
    optr = buffer2;
    
    // execute remaining push stages
    qsm::seq_cgrqueue<GI, VT> inq;
    for(UGI s = 1; s < npushstages; ++s) {
      UGI const nparts = schedule.get_num_push_parts(s);
      for(UGI p = 0; p < nparts; ++p) {
        ispace<UGI> const & part_space = schedule.get_push_part_space(s, p);
        ispace<UGI> const & subparts = schedule.get_subpart_range(s, p);
        
        VT const * inptr = &iptr[part_space.start];
        inq.configure(inptr);
        
        for(UGI sp = subparts.start, idx = 0; sp < subparts.end; ++sp, ++idx) {
          ispace<UGI> const & subpart_space = schedule.get_push_part_space(s+1, sp);
          outq[idx].configure(&optr[subpart_space.start]);
        }
        
        DI const * push_map = schedule.get_push_map(s, p);
        
        for(UGI i = 0; i < part_space.size; i+=BufferSize) {
          #pragma novector
          for(UGI k = 0; k < subparts.size; ++k) {
            bufptrs[k] = outq[k].get();
            bufidx[k] = 0;
          }
          
          int niter = std::min((int)(part_space.size-i), BufferSize);
          
          #pragma loop_count max(BufferSize)
          for(int j = 0; j < niter; ++j) {
            auto val = inq.pop();
            DI nbkts = *push_map++;
            for(DI k = 0; k < nbkts; ++k) {
              bufptrs[*push_map][bufidx[*push_map]] = val; //inptr[i+j];
              bufidx[*push_map]++;
              push_map++;
            }
          }
          
          for(UGI k = 0; k < subparts.size; ++k) {
            #pragma forceinline recursive
            outq[k].flush(bufidx[k]);
          }
        }
        
        for(UGI k = 0; k < subparts.size; ++k) {
          outq[k].deactivate();
        }
      }
      
      VT * tptr = const_cast<VT *>(iptr);
      iptr = optr;
      optr = tptr;
      
      _mm_mfence();
    }
  } else {
    optr = buffer1;
    
    qsm::seq_cgrqueue<GI, VT> inq;
    
    qsm::cgwqueue<GI, VT, BufferSize> outq;
    VT * bufptr;
    
    outq.configure(optr);
    
    for(GI i = 0; i < ninspc; ++i) {
      inq.configure(inspc[i].start, input);
      
      for(GI j = 0; j < inspc[i].size; j+=BufferSize) {
        bufptr = outq.get();
        int niter = std::min(inspc[i].size-j, BufferSize);
        for(int l = 0; l < niter; ++l) {
          bufptr[l] = inq.pop();//gather(inq);
          //inq.advance();
        }
        
        outq.flush(niter);
      }
      
      outq.deactivate();
    }
    
    iptr = optr;
    
    _mm_mfence();
  }
  
  // pull phase
  {
    UGI const npullparts = schedule.get_num_pull_parts();
    
    qsm::cgwqueue_pack<BufferSize, GI, T...> outq;
    std::tuple<T *...> bufptr;
    
    outq.configure(output);
    
    for(UGI p = 0; p < npullparts; ++p) {
      ispace<UGI> const & inspace = schedule.get_push_part_space(npushstages, p);
      ispace<UGI> const & outspace = schedule.get_pull_part_space(p);
      
      UGI const * pull_map = schedule.get_pull_map(p);
      VT const * inptr = &iptr[inspace.start];
      
      for(UGI i = 0; i < outspace.size; i+=BufferSize) {
        bufptr = outq.get();
        
        int niter = std::min((int)(outspace.size-i), BufferSize);
        
        #pragma loop_count max(BufferSize)
        for(int j = 0; j < niter; ++j) {
          #pragma forceinline
          scatter<GI, T...>(bufptr, j, inptr[pull_map[i+j]]);
        }
        
        outq.flush(niter);
      }
    }
    
    outq.deactivate();
  }
}

template<typename GI, typename DI>
struct distribution_schedule_helper {
public:
  /*static void get_entity_and_nbuckets_count(
    int const nparts,
    distribution_schedule<GI, DI> const * schedules,
    GI * entity_count, GI * nbkts_count
  ) {
    using UGI = typename distribution_schedule<GI, DI>::uindex_type;
    
    *entity_count = 0;
    *nbkts_count = 0;
    
    for(int t = 0; t < nparts; ++t) {
      GI const ninspc = schedules[t].get_num_input_spaces();
      ispace<GI> const * inspc = schedules[t].get_input_spaces();
      
      UGI const npushstages = schedules[t].get_num_push_stages();
      if(npushstages > 0) {
        {
          DI const * push_map = schedules[t].get_push_map(0, 0);
          for(GI i = 0; i < ninspc; ++i) {
            *entity_count += inspc[i].size;
            for(GI j = 0; j < inspc[i].size; ++j) {
              DI nbkts = *push_map++;
              *nbkts_count += nbkts;
              push_map += nbkts;
            }
          }
        }
        
        for(UGI s = 1; s < npushstages; ++s) {
          UGI const nparts = schedules[t].get_num_push_parts(s);
          for(UGI p = 0; p < nparts; ++p) {
            ispace<UGI> const & part_space = schedules[t].get_push_part_space(s, p);
            *entity_count += part_space.size;
            DI const * push_map = schedules[t].get_push_map(s, p);
            for(UGI i = 0; i < part_space.size; ++i) {
              DI nbkts = *push_map++;
              *nbkts_count += nbkts;
              push_map += nbkts;
            }
          }
        }
      }
    }
  }*/
  
  template<typename... T>
  static
  size_t get_distribution_memsize(
    int const nparts,
    distribution_schedule<GI, DI> const * schedules,
    size_t * insizestore = nullptr, size_t * outsizestore = nullptr,
    size_t * mapsizestore = nullptr
  ) {
    using VT = std::tuple<T...>;
    using UGI = typename distribution_schedule<GI, DI>::uindex_type;
    
    size_t sz1 = qsm::sum_sizes(VT());
    size_t sz2 = sizeof(VT);
    
    size_t outsize = 0;
    size_t insize = 0;
    size_t mapsize = 0;
    
    for(int p = 0; p < nparts; ++p) {
      GI const ninspc = schedules[p].get_num_input_spaces();
      ispace<GI> const * inspc = schedules[p].get_input_spaces();
      
      UGI const npushstages = schedules[p].get_num_push_stages();
      
      // push phase
      if(npushstages > 0) {
        // first push stage
        {
          for(GI i = 0; i < ninspc; ++i) {
            insize += sz1*inspc[i].size;
          }
          
          ispace<UGI> const & subparts = schedules[p].get_subpart_range(0, 0);
          for(UGI sp = subparts.start; sp < subparts.end; ++sp) {
            ispace<UGI> const & subpart_space = schedules[p].get_push_part_space(1, sp);
            outsize += sz2*subpart_space.size;
          }
          
          mapsize += sizeof(DI)*schedules[p].get_push_map_size(0, 0);
        }
        
        // remaining push stages
        for(UGI s = 1; s < npushstages; ++s) {
          UGI const nparts = schedules[p].get_num_push_parts(s);
          for(UGI pp = 0; pp < nparts; ++pp) {
            ispace<UGI> const & part_space = schedules[p].get_push_part_space(s, pp);
            insize += sz2*part_space.size;
            
            ispace<UGI> const & subpart_range = schedules[p].get_subpart_range(s, pp);
            for(UGI sp = subpart_range.start; sp < subpart_range.end; ++sp) {
              ispace<UGI> const & subpart_space = schedules[p].get_push_part_space(s+1, sp);
              outsize += sz2*subpart_space.size;
            }
            
            mapsize += sizeof(DI)*schedules[p].get_push_map_size(s, pp);
          }
        }
      } else {
        for(GI i = 0, idx = 0; i < ninspc; ++i) {
          insize += sz1*inspc[i].size;
          outsize += sz2*inspc[i].size;
        }
      }
      
      // pull phase
      UGI const npullparts = schedules[p].get_num_pull_parts();
      for(UGI pp = 0; pp < npullparts; ++pp) {
        ispace<UGI> const & part_space = schedules[p].get_pull_part_space(pp);
        insize += sz2*part_space.size;
        outsize += sz1*part_space.size;
        mapsize += sizeof(UGI)*schedules[p].get_pull_map_size(pp);
      }
    }
    
    if(insizestore != nullptr)
      *insizestore = insize;
    if(outsizestore != nullptr)
      *outsizestore = outsize;
    if(mapsizestore != nullptr)
      *mapsizestore = mapsize;
    
    return insize+outsize+mapsize;
  }
  
  template<typename... T>
  static
  size_t get_distribution_memsize3(
    int const nparts,
    distribution_schedule<GI, DI> const * schedules,
    size_t * insizestore = nullptr, size_t * outsizestore = nullptr,
    size_t * mapsizestore = nullptr
  ) {
    using VT = std::tuple<T...>;
    using UGI = typename distribution_schedule<GI, DI>::uindex_type;
    
    size_t sz1 = qsm::sum_sizes(VT());
    size_t sz2 = sizeof(VT);
    
    size_t outsize = 0;
    size_t insize = 0;
    size_t mapsize = 0;
    
    for(int p = 0; p < nparts; ++p) {
      GI const ninspc = schedules[p].get_num_input_spaces();
      ispace<GI> const * inspc = schedules[p].get_input_spaces();
      
      UGI const npushstages = schedules[p].get_num_push_stages();
      
      // push phase
      if(npushstages > 0) {
        // first push stage
        {
          for(GI i = 0; i < ninspc; ++i) {
            insize += sz2*inspc[i].size;
          }
          
          ispace<UGI> const & subparts = schedules[p].get_subpart_range(0, 0);
          for(UGI sp = subparts.start; sp < subparts.end; ++sp) {
            ispace<UGI> const & subpart_space = schedules[p].get_push_part_space(1, sp);
            outsize += sz2*subpart_space.size;
          }
          
          mapsize += sizeof(DI)*schedules[p].get_push_map_size(0, 0);
        }
        
        // remaining push stages
        for(UGI s = 1; s < npushstages; ++s) {
          UGI const nparts = schedules[p].get_num_push_parts(s);
          for(UGI pp = 0; pp < nparts; ++pp) {
            ispace<UGI> const & part_space = schedules[p].get_push_part_space(s, pp);
            insize += sz2*part_space.size;
            
            ispace<UGI> const & subpart_range = schedules[p].get_subpart_range(s, pp);
            for(UGI sp = subpart_range.start; sp < subpart_range.end; ++sp) {
              ispace<UGI> const & subpart_space = schedules[p].get_push_part_space(s+1, sp);
              outsize += sz2*subpart_space.size;
            }
            
            mapsize += sizeof(DI)*schedules[p].get_push_map_size(s, pp);
          }
        }
      } else {
        for(GI i = 0, idx = 0; i < ninspc; ++i) {
          insize += sz2*inspc[i].size;
          outsize += sz2*inspc[i].size;
        }
      }
      
      // pull phase
      UGI const npullparts = schedules[p].get_num_pull_parts();
      for(UGI pp = 0; pp < npullparts; ++pp) {
        ispace<UGI> const & part_space = schedules[p].get_pull_part_space(pp);
        insize += sz2*part_space.size;
        outsize += sz1*part_space.size;
        mapsize += sizeof(UGI)*schedules[p].get_pull_map_size(pp);
      }
    }
    
    if(insizestore != nullptr)
      *insizestore = insize;
    if(outsizestore != nullptr)
      *outsizestore = outsize;
    if(mapsizestore != nullptr)
      *mapsizestore = mapsize;
    
    return insize+outsize+mapsize;
  }
  
  template<typename T>
  static
  size_t get_distribution_memsize(
    int const nparts,
    distribution_schedule<GI, DI> const * schedules,
    size_t * insizestore = nullptr, size_t * outsizestore = nullptr,
    size_t * mapsizestore = nullptr
  ) {
    using UGI = typename distribution_schedule<GI, DI>::uindex_type;
    
    size_t outsize = 0;
    size_t insize = 0;
    size_t mapsize = 0;
    
    for(int p = 0; p < nparts; ++p) {
      GI const ninspc = schedules[p].get_num_input_spaces();
      ispace<GI> const * inspc = schedules[p].get_input_spaces();
      
      UGI const npushstages = schedules[p].get_num_push_stages();
      
      // push phase
      if(npushstages > 0) {
        // first push stage
        {
          for(GI i = 0; i < ninspc; ++i) {
            insize += sizeof(T)*inspc[i].size;
          }
          
          ispace<UGI> const & subparts = schedules[p].get_subpart_range(0, 0);
          for(UGI sp = subparts.start; sp < subparts.end; ++sp) {
            ispace<UGI> const & subpart_space = schedules[p].get_push_part_space(1, sp);
            outsize += sizeof(T)*subpart_space.size;
          }
          
          mapsize += sizeof(DI)*schedules[p].get_push_map_size(0, 0);
        }
        
        // remaining push stages
        for(UGI s = 1; s < npushstages; ++s) {
          UGI const nparts = schedules[p].get_num_push_parts(s);
          for(UGI pp = 0; pp < nparts; ++pp) {
            ispace<UGI> const & part_space = schedules[p].get_push_part_space(s, pp);
            insize += sizeof(T)*part_space.size;
            
            ispace<UGI> const & subpart_range = schedules[p].get_subpart_range(s, pp);
            for(UGI sp = subpart_range.start; sp < subpart_range.end; ++sp) {
              ispace<UGI> const & subpart_space = schedules[p].get_push_part_space(s+1, sp);
              outsize += sizeof(T)*subpart_space.size;
            }
            
            mapsize += sizeof(DI)*schedules[p].get_push_map_size(s, pp);
          }
        }
      } else {
        for(GI i = 0, idx = 0; i < ninspc; ++i) {
          insize += sizeof(T)*inspc[i].size;
          outsize += sizeof(T)*inspc[i].size;
        }
      }
      
      // pull phase
      UGI const npullparts = schedules[p].get_num_pull_parts();
      for(UGI pp = 0; pp < npullparts; ++pp) {
        ispace<UGI> const & part_space = schedules[p].get_pull_part_space(pp);
        outsize += sizeof(T)*part_space.size;
        insize += sizeof(T)*part_space.size;
        mapsize += sizeof(UGI)*schedules[p].get_pull_map_size(pp);
      }
    }
    
    if(insizestore != nullptr)
      *insizestore = insize;
    if(outsizestore != nullptr)
      *outsizestore = outsize;
    if(mapsizestore != nullptr)
      *mapsizestore = mapsize;
    
    return insize+outsize+mapsize;
  }
};

