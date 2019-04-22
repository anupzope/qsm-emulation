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
#include "spmatrix.hpp"
#include "types.hpp"

template<typename MT, typename VT>
struct spmv_traits {
  static size_t get_flop_count_per_nnz() {
    static_assert(false, "Not implemented");
  }
};

template<>
struct spmv_traits<double, double> {
  static size_t get_flop_count_per_nnz() {
    return 2;
  }
};

template<typename T, int R, int CR, int C>
struct spmv_traits<matrix<T, R, CR>, matrix<T, CR, C>> {
  static size_t get_flop_count_per_nnz() {
    return R*CR*C*2 + R*C;
  }
};

// Usage:
// int retval = 0;
// init()
// if((retval = preprocess()) == 0) {
//   execute()
//   calculate_metrics()
//   postprocess()
// } else {
//   //print error
// }
// finalize()

template<typename GI, typename MT, typename VT>
class spmv_test {
protected:
  // input parameters
  GI nrepeat;
  GI nparts;
  ispace<GI> const * parts;
  spmatrixstructure_crs<GI> const * matstr;
  MT const * a;
  VT const * x0;
  VT * x1;
  
  // internal data structures
  std::vector<double> time;
  std::string errmsg;
  
  // metrics
  size_t nflops;
  size_t data_volume;
  double mintime;
  double maxtime;
  double avgtime;
  double flop2byte;
  double bandwidth_mbps;
  
protected:
template<size_t ... I>
  size_t calculate_flop_count_per_nnz() const {
    return spmv_traits<MT, VT>::get_flop_count_per_nnz();
  }
  
public:
  spmv_test() : matstr(nullptr), a(nullptr), x0(nullptr), x1(nullptr),
  nrepeat(0), nflops(0), data_volume(0), mintime(0.0), maxtime(0.0),
  avgtime(0.0), flop2byte(0.0), bandwidth_mbps(0.0) {
  }
  
  virtual ~spmv_test() {}
  
  std::string get_error_message() const {
    return errmsg;
  }
  
  size_t get_flops_count() const {
    return nflops;
  }
  
  size_t get_data_volume() const {
    return data_volume;
  }
  
  double get_min_time() const {
   return mintime;
  }
  
  double get_max_time() const {
   return maxtime;
  }
  
  double get_avg_time() const {
   return avgtime;
  }
  
  double get_flop2byte() const {
    return flop2byte;
  }
  
  double get_bandwidth() const {
    return bandwidth_mbps;
  }
  
  void init(
    int const nrepeat, GI const nparts, ispace<GI> const * parts,
    spmatrixstructure_crs<GI> const * matstr,
    VT * x1, MT const * a, VT const * x0
  ) {
    this->nrepeat = nrepeat;
    this->nparts = nparts;
    this->parts = parts;
    this->matstr = matstr;
    this->x1 = x1;
    this->a = a;
    this->x0 = x0;
    
    time.resize(nrepeat);
    for(int r = 0; r < nrepeat; ++r)
      time[r] = 0.0;
    
    nflops = 0;
    data_volume = 0;
    mintime = 0.0;
    maxtime = 0.0;
    avgtime = 0.0;
    flop2byte = 0.0;
    bandwidth_mbps = 0.0;
  }
  
  virtual void finalize() {}
  
  virtual int preprocess() {
    return 0;
  }
  
  virtual void postprocess() {}
  
  virtual void execute() {
    timespec start_time, end_time;
    
    for(int r = 0; r < nrepeat; ++r) {
      clock_gettime(CLOCK_MONOTONIC, &start_time);
      this->execute_internal();
      clock_gettime(CLOCK_MONOTONIC, &end_time);
      
      time[r] = double(end_time.tv_sec-start_time.tv_sec) +
        1e-9*(end_time.tv_nsec-start_time.tv_nsec);
    }
  }
  
  virtual void calculate_metrics() {
    this->calculate_preliminary_metrics();
    
    mintime = std::numeric_limits<double>::max();
    maxtime = 0.0;
    avgtime = 0.0;
    for(int r = 1; r < nrepeat; ++r) {
      mintime = std::min(mintime, time[r]);
      maxtime = std::max(maxtime, time[r]);
      avgtime += time[r];
    }
    avgtime /= (double)(nrepeat-1);
    
    double data_volume_mb = (double)data_volume/1024.0/1024.0;
    bandwidth_mbps = data_volume_mb/mintime;
    flop2byte = (double)nflops/(double)data_volume;
  }
  
  virtual void execute_internal() = 0;
  virtual void calculate_preliminary_metrics() = 0;
};

class spmv_pack_helper {
  template<typename RT, typename ... LT, size_t ... I>
  static
  void assign1(
    std::tuple<LT...> & lhs, RT const & rhs, std::index_sequence<I...>
  ) {
    using swallow = int[];
    (void)swallow{0, (void(std::get<I>(lhs) = rhs), 0) ... };
  }
  
  template<typename GI, typename VT, typename ... PT, size_t ... I>
  static void assign1(
    std::tuple<PT * ...> const & ptrs, GI const idx,
    VT const & val, std::index_sequence<I...>
  ) {
    using swallow = int[];
    (void)swallow{0, (void(std::get<I>(ptrs)[idx] = val), 0) ... };
  }
  
  template<typename GI, typename ... PT, size_t ... I>
  static void assign1(
    std::tuple<PT * ...> const & ptrs, GI const idx,
    std::tuple<PT ...> const & val, std::index_sequence<I...>
  ) {
    using swallow = int[];
    (void)swallow{0, (void(std::get<I>(ptrs)[idx] = std::get<I>(val)), 0) ... };
  }
  
  template<typename GI, typename ... PT, size_t ... I>
  static
  std::ostream & print1(
    std::ostream & stream, std::tuple<PT * ...> const & ptrs, GI const idx,
    std::index_sequence<I...>
  ) {
    using swallow = int[];
    (void)swallow{0, (void(stream << (I == 0 ? "" : " ") << std::get<I>(ptrs)[idx]), 0)... };
    return stream;
  }
  
  template<typename GI, typename MT, typename ... VT, size_t ... I>
  static
  void multiply_add1(
    std::tuple<VT ...> & sum,
    MT const & a,
    std::tuple<VT const * ...> x0,
    GI idx,
    std::index_sequence<I...>
  ) {
    auto r = { (std::get<I>(sum) += a * std::get<I>(x0)[idx]) ... };
  }
  
public:
  template<typename RT, typename ... LT>
  static
  void assign(
    std::tuple<LT...> & lhs, RT const & rhs
  ) {
    #pragma forceinline
    assign1(lhs, rhs, std::make_index_sequence<sizeof...(LT)>{});
  }
  
  template<typename GI, typename VT, typename ... PT>
  static
  void assign(
    std::tuple<PT * ...> const & ptrs, GI const idx,
    VT const & val
  ) {
    #pragma forceinline
    assign1(ptrs, idx, val, std::make_index_sequence<sizeof...(PT)>{});
  }
  
  template<typename GI, typename ... PT>
  static
  void assign(
    std::tuple<PT * ...> const & ptrs, GI const idx,
    std::tuple<PT ...> const & val
  ) {
    #pragma forceinline
    assign1(ptrs, idx, val, std::make_index_sequence<sizeof...(PT)>{});
  }
  
  template<typename GI, typename ... PT>
  static
  std::ostream & print(
    std::ostream & stream, std::tuple<PT * ...> const & ptrs, GI const idx
  ) {
    #pragma forceinline
    return print1(stream, ptrs, idx, std::make_index_sequence<sizeof...(PT)>{});
  }

  template<typename GI, typename MT, typename ... VT>
  static
  void multiply_add( // sum += a * x0[idx]
    std::tuple<VT ...> & sum,
    MT const & a,
    std::tuple<VT const * ...> x0,
    GI idx
  ) {
    #pragma forceinline
    multiply_add1(
      sum, a, x0, idx, std::make_index_sequence<sizeof...(VT)>{}
    );
  }
};

template<typename GI, typename MT, typename ... VT>
class spmv_pack_test {
protected:
  // input parameters
  GI nrepeat;
  GI nparts;
  ispace<GI> const * parts;
  spmatrixstructure_crs<GI> const * matstr;
  MT const * a;
  std::tuple<VT const * ...> x0;
  std::tuple<VT * ...> x1;
  
  // internal data structures
  std::vector<double> time;
  std::string errmsg;
  
  // metrics
  size_t nflops;
  size_t data_volume;
  double mintime;
  double maxtime;
  double avgtime;
  double flop2byte;
  double bandwidth_mbps;
  
protected:
  template<size_t ... I>
  size_t calculate_flop_count_per_nnz(std::index_sequence<I...>) const {
    size_t flop_count = 0;
    auto r = {
      (
        flop_count +=
          spmv_traits<
            MT,
            typename std::tuple_element<I, std::tuple<VT...>>::type
          >::get_flop_count_per_nnz()
      )...
    };
    return flop_count;
  }
  
public:
  spmv_pack_test() : matstr(nullptr), a(nullptr), x0{}, x1{},
  nrepeat(0), nflops(0), data_volume(0), mintime(0.0), maxtime(0.0),
  avgtime(0.0), flop2byte(0.0), bandwidth_mbps(0.0) {
  }
  
  virtual ~spmv_pack_test() {}
  
  std::string get_error_message() const {
    return errmsg;
  }
  
  size_t get_flops_count() const {
    return nflops;
  }
  
  size_t get_data_volume() const {
    return data_volume;
  }
  
  double get_min_time() const {
   return mintime;
  }
  
  double get_max_time() const {
   return maxtime;
  }
  
  double get_avg_time() const {
   return avgtime;
  }
  
  double get_flop2byte() const {
    return flop2byte;
  }
  
  double get_bandwidth() const {
    return bandwidth_mbps;
  }
  
  void init(
    int const nrepeat, GI const nparts, ispace<GI> const * parts,
    spmatrixstructure_crs<GI> const * matstr,
    std::tuple<VT * ...> const & x1, MT const * a,
    std::tuple<VT const * ...> const & x0
  ) {
    this->nrepeat = nrepeat;
    this->nparts = nparts;
    this->parts = parts;
    this->matstr = matstr;
    this->x1 = x1;
    this->a = a;
    this->x0 = x0;
    
    time.resize(nrepeat);
    for(int r = 0; r < nrepeat; ++r)
      time[r] = 0.0;
    
    nflops = 0;
    data_volume = 0;
    mintime = 0.0;
    maxtime = 0.0;
    avgtime = 0.0;
    flop2byte = 0.0;
    bandwidth_mbps = 0.0;
  }
  
  virtual void finalize() {}
  
  virtual int preprocess() {
    return 0;
  }
  
  virtual void postprocess() {}
  
  virtual void execute() {
    timespec start_time, end_time;
    
    for(int r = 0; r < nrepeat; ++r) {
      clock_gettime(CLOCK_MONOTONIC, &start_time);
      this->execute_internal();
      clock_gettime(CLOCK_MONOTONIC, &end_time);
      
      time[r] = double(end_time.tv_sec-start_time.tv_sec) +
        1e-9*(end_time.tv_nsec-start_time.tv_nsec);
    }
  }
  
  virtual void calculate_metrics() {
    this->calculate_preliminary_metrics();
    
    mintime = std::numeric_limits<double>::max();
    maxtime = 0.0;
    avgtime = 0.0;
    for(int r = 1; r < nrepeat; ++r) {
      mintime = std::min(mintime, time[r]);
      maxtime = std::max(maxtime, time[r]);
      avgtime += time[r];
    }
    avgtime /= (double)(nrepeat-1);
    
    double data_volume_mb = (double)data_volume/1024.0/1024.0;
    bandwidth_mbps = data_volume_mb/mintime;
    flop2byte = (double)nflops/(double)data_volume;
  }
  
  virtual void execute_internal() = 0;
  virtual void calculate_preliminary_metrics() = 0;
};
