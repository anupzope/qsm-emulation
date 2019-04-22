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

#include "spmatrix.hpp"
#include "spmv_test_temporal_v1.hpp"
#include "spmv_test_qsm0_v1.hpp"
#include "spmv_test_qsm0_v2.hpp"
#include "types.hpp"

#include <cstdio>
#include <random>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <cstring>

#include <argp.h>
#include <fenv.h>
#include <signal.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef USING_ITT
#include <ittnotify.h>
#endif

template<typename GI>
void simple_partition(
  GI const nparts,
  GI const nrows,
  GI const * rowoffs,
  GI * parts
) {
  GI qut = nrows/nparts;
  GI rem = nrows%nparts;
  for(GI p = 0; p < nparts; ++p)
    parts[p+1] = parts[p] + qut + (p < rem ? 1 : 0);
}

template<typename GI>
void outdegree_weighted_partition(
  GI const nrows,
  GI const * rowoffs,
  GI const nparts,
  GI * partoffs
) {
  GI const nnz = rowoffs[nrows];
  GI qut = nnz/nparts;
  GI rem = nnz%nparts;
  
  partoffs[0] = 0;
  for(GI p = 0; p < nparts; ++p)
    partoffs[p+1] = partoffs[p] + qut + (p < rem ? 1 : 0);
  
  for(GI r = 0, p = 1; r < nrows && p < nparts; ++r) {
    if(rowoffs[r+1] >= partoffs[p]) {
      partoffs[p] = r+1;
      ++p;
    }
  }
  partoffs[nparts] = nrows;
}

void fpe_handler(int sig) {
  fprintf(stderr, "Floating point exception.\n");
  exit(1);
}

enum version {
  SPMV_TEMPORAL,
  SPMV_QSM0_V1,
  SPMV_QSM0_V2
};

static char const * version_names[] = {
  "temporal",
  "qsm0v1",
  "qsm0v2"
};

template<typename GI, typename LI>
struct program_options {
  char const * casename;
  GI version;
  bool out;
  GI lookback;
  GI lookahead;
  typename std::make_unsigned<GI>::type ntailbits = 12;
  typename std::make_unsigned<GI>::type nbitsperstage = 6;
  GI nrepeat;
  char const * batchfilename;
  
  program_options() {
    casename = nullptr;
    version = SPMV_TEMPORAL;
    out = false;
    lookback = 8192;
    lookahead = 4096;
    ntailbits = 12;
    nbitsperstage = 6;
    nrepeat = 11;
    batchfilename = nullptr;
  }
};

template<typename GI, typename LI>
static int parse_options(int key, char * arg, argp_state * state) {
  program_options<GI, LI> * opt = (program_options<GI, LI>*)state->input;
  switch(key) {
  case 401:
    {
      char * end;
      errno = 0;
      unsigned long version = std::strtoul(arg, &end, 10);
      if(errno == ERANGE) {
        errno = 0;
        argp_failure(state, 1, 0, "version is out of range");
      }
      if(version > 2) {
        argp_failure(state, 1, 0, "version must be between [0, 2]");
      }
      opt->version = version;
    }
    break;
  case 402:
    {
      opt->out = true;
    }
    break;
  case 403:
    {
      char * end;
      errno = 0;
      unsigned long lookback = std::strtoul(arg, &end, 10);
      if(errno == ERANGE) {
        errno = 0;
        argp_failure(state, 1, 0, "lookback is out of range");
      }
      if(lookback >= std::numeric_limits<LI>::max()) {
        argp_failure(state, 1, 0, "lookback is out of range");
      }
      opt->lookback = lookback;
    }
    break;
  case 404:
    {
      char * end;
      errno = 0;
      unsigned long lookahead = std::strtoul(arg, &end, 10);
      if(errno == ERANGE) {
        errno = 0;
        argp_failure(state, 1, 0, "lookahead is out of range");
      }
      if(lookahead >= std::numeric_limits<LI>::max()) {
        argp_failure(state, 1, 0, "lookahead is out of range");
      }
      opt->lookahead = lookahead;
    }
    break;
  case 405:
    {
      char * end;
      errno = 0;
      unsigned long value = std::strtoul(arg, &end, 10);
      if(errno == ERANGE) {
        errno = 0;
        argp_failure(state, 1, 0, "nbitsperstage is out of range");
      }
      if(value > 64) {
        argp_failure(state, 1, 0, "nbitsperstage is out of range");
      }
      
      opt->nbitsperstage = value;
    }
    break;
  case 406:
    {
      char * end;
      errno = 0;
      unsigned long value = std::strtoul(arg, &end, 10);
      if(errno == ERANGE) {
        errno = 0;
        argp_failure(state, 1, 0, "ntailbits is out of range");
      }
      if(value > 64) {
        argp_failure(state, 1, 0, "ntailbits is out of range");
      }
      
      opt->ntailbits = value;
    }
    break;
  case 407:
    {
      char * end;
      errno = 0;
      unsigned long value = std::strtoul(arg, &end, 10);
      if(errno == ERANGE) {
        errno = 0;
        argp_failure(state, 1, 0, "nrepeat is out of range");
      }
      if(value < 2 || value >= 256) {
        argp_failure(state, 1, 0, "nrepeat must be between [2, 256)");
      }
      opt->nrepeat = value;
    }
    break;
  case 408:
    {
      if(std::strlen(arg) > 0)
        opt->batchfilename = arg;
    }
    break;
  case ARGP_KEY_ARG:
    {
      switch(state->arg_num) {
      case 0:
        opt->casename = arg;
        break;
      default:
        // too many arguments
        argp_usage(state);
      }
    }
    break;
  case ARGP_KEY_END:
    {
      if(state->arg_num < 1)
        // not enough arguments
        argp_usage(state);
    }
    break;
  default:
    return ARGP_ERR_UNKNOWN;
  }
  return 0;
}

template<typename ... VT>
struct pack_traits {
  using store_type = pack_store<VT...>;
  
  using ptrs_type = std::tuple<VT * ...>;
  
  using constptrs_type = std::tuple<VT const * ...>;
  
  template<typename GI, typename MT>
  using test_base_type = spmv_pack_test<GI, MT, VT...>;
  
  template<typename GI, typename MT>
  using test_temporal_v1_type = spmv_pack_test_temporal_v1<GI, MT, VT...>;
  
  template<typename GI, typename LI, typename DI, typename MT>
  using test_qsm0_v2_type = spmv_pack_test_qsm0_v2<GI, LI, DI, MT, VT...>;
};

int main(int argc, char * argv[]) {
#ifdef _OPENMP
  omp_set_dynamic(0);
  omp_set_nested(0);
  int const maxthreads = omp_get_max_threads();
  int numthreads = 0;
  #pragma omp parallel
  {
    #pragma omp atomic update
    numthreads += 1;
  }
  if(maxthreads != numthreads) {
    printf("The number of threads counted not the same as the number of threads requested\n");
    return 1;
  }
  char const * env_omp_num_threads = std::getenv("OMP_NUM_THREADS");
  char const * env_kmp_affinity = std::getenv("KMP_AFFINITY");
#else
  int const maxthreads = 1;
  int numthreads = 1;
#endif
  
  using GI = int;
  using LI = short;
  using DI = uint8_t;
  using VT = matrix<double, 5, 1>;
  using MT = matrix<double, 5, 5>;
  
  feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW  | FE_UNDERFLOW);
  signal(SIGFPE, fpe_handler);
  
  program_options<GI, LI> opt;
  
  struct argp_option options[] = {
    {"version", 401, "NUM", 0, "NUM=0: temporal, NUM=1: qsm0v1, NUM=2: qsm0v2"},
    {"out", 402, 0, 0, "print the result or not"},
    {"lookback", 403, "NUM", 0, "look-back distance"},
    {"lookahead", 404, "NUM", 0, "look-ahead distance"},
    {"nbitsperstage", 405, "NUM", 0, "number of bits processed per stage of the distribution schedule"},
    {"ntailbits", 406, "NUM", 0, "number of trailing bits in the distribution schedule"},
    {"nrepeat", 407, "NUM", 0, "number of times to repeat the experiment"},
    {"batch", 408, "STRING", 0, "file name for appending results"},
    {0}
  };
  
  argp ap = {options, parse_options<GI, LI>, "CASE-NAME"};
  
  int parserv = 0;
  if((parserv = argp_parse(&ap, argc, argv, 0, 0, &opt)))
    return parserv;
  
  std::string solfilename(opt.casename);
  solfilename += "-";
  solfilename += version_names[opt.version];
  solfilename += ".sol";
  
  printf("casename          : %s\n", opt.casename);
  printf("version           : %s\n", version_names[opt.version]);
  printf("lookback          : %d\n", opt.lookback);
  printf("lookahead         : %d\n", opt.lookahead);
  printf("nbitsperstage     : %d\n", opt.nbitsperstage);
  printf("ntailbits         : %d\n", opt.ntailbits);
  printf("nrepeat           : %d\n", opt.nrepeat);
  if(opt.out)
  printf("output file       : %s\n", solfilename.c_str());
  printf("number of threads : %d\n", numthreads);
  printf("OMP_NUM_THREADS   : %s\n", env_omp_num_threads);
  printf("KMP_AFFINITY      : %s\n", env_kmp_affinity);
  
  std::string crsfilename(opt.casename);
  crsfilename += ".crs";
  
  spmatrixstructure_crs<GI> * matstr = new spmatrixstructure_crs<GI>();
  
  std::ifstream crsfile(crsfilename.c_str());
  if(!crsfile) {
    fprintf(stderr, "could not open file %s\n", crsfilename.c_str());
    return 1;
  }
  crsfile >> matstr->nrows >> matstr->ncols >> matstr->nnz;
  printf("Matrix info       : nrows = %d, ncols = %d, nnz = %d\n",
    matstr->nrows, matstr->ncols, matstr->nnz);
  if(!crsfile) {
    fprintf(stderr, "could not read matrix metadata\n");
    return 1;
  }
  
  matstr->rowoffs = new GI[matstr->nrows+1];
  matstr->colidxs = new GI[matstr->nnz];
  matstr->rowoffs[0] = 0;
  for(GI r = 0, idx = 0; r < matstr->nrows; ++r) {
    GI noffs;
    crsfile >> noffs;
    if(!crsfile) {
      fprintf(stderr, "error reading number of row offsets at %s:%d\n",
        crsfilename.c_str(), r+1);
      crsfile.close();
      return 1;
    }
    for(GI i = 0; i < noffs; ++i) {
      crsfile >> matstr->colidxs[idx++];
      if(!crsfile) {
        fprintf(stderr, "error reading column index at %d th offset at %s:%d\n",
          i, crsfilename.c_str(), r+1);
        crsfile.close();
        return 1;
      }
    }
    matstr->rowoffs[r+1] = matstr->rowoffs[r] + noffs;
  }
  crsfile.close();
  
  std::vector<GI> partoffs(numthreads+1);
  printf("partitioning index space\n");
  outdegree_weighted_partition(
    matstr->nrows, matstr->rowoffs, numthreads, &partoffs[0]
  );
  
  std::vector<ispace<GI>> parts(numthreads);
  for(GI p = 0; p < numthreads; ++p) { 
    parts[p] = ispace<GI>(partoffs[p], partoffs[p+1]);
  }
  
  size_t nflops = 0;
  size_t data_volume = 0;
  double mintime = 0;
  double maxtime = 0;
  double avgtime = 0;
  double flop2byte = 0;
  double bandwidth_mbps = 0;
  bool test_complete = false;
  
#define SPMV_PACK_TEST
#ifdef SPMV_PACK_TEST
  using PT = pack_traits<VT, VT, VT, VT, VT>;
  PT::store_type x0, x1;
  std::vector<MT> a(matstr->nnz);
  
  x1.allocate(matstr->nrows);
  x0.allocate(matstr->ncols);
  PT::ptrs_type x1ptrs = x1.get();
  PT::ptrs_type x0ptrs = x0.get();
  PT::constptrs_type x0cptrs = x0.get_const();
  
  for(GI r = 0; r < matstr->nrows; ++r) {
    spmv_pack_helper::assign(x1ptrs, r, 0.0);
  }
  std::mt19937 rng;
  std::uniform_real_distribution<> dis(1.0, 10.0);
  for(GI c = 0; c < matstr->ncols; ++c) {
    spmv_pack_helper::assign(x0ptrs, c, dis(rng));
  }
  for(GI i = 0; i < matstr->nnz; ++i) {
    a[i] = 1.0;
  }
  
  PT::template test_base_type<GI, MT> * test = nullptr;
  switch(opt.version) {
  case SPMV_TEMPORAL:
    test = new PT::template test_temporal_v1_type<GI, MT>();
    break;
  case SPMV_QSM0_V1:
    {
      printf("Pack version of qsm0v1 not implemented.\n");
      exit(1);
    }
    break;
  case SPMV_QSM0_V2:
    {
      auto obj = new PT::template test_qsm0_v2_type<GI, LI, DI, MT>();
      obj->set_parameters(
        opt.lookahead, opt.lookback, opt.ntailbits, opt.nbitsperstage
      );
      test = obj;
    }
    break;
  }
  
  test->init(
    opt.nrepeat, numthreads, &parts[0],
    matstr, x1ptrs, &a[0], x0cptrs
  );
  int retval = test->preprocess();
  if(retval == 0) {
#ifdef USING_ITT
    __itt_resume();
#endif
    test->execute();
#ifdef USING_ITT
    __itt_pause();
#endif
    
    test->calculate_metrics();
    
    nflops = test->get_flops_count();
    data_volume = test->get_data_volume();
    mintime = test->get_min_time();
    maxtime = test->get_max_time();
    avgtime = test->get_avg_time();
    flop2byte = test->get_flop2byte();
    bandwidth_mbps = test->get_bandwidth();
    test_complete = true;
    
    test->postprocess();
  } else {
    fprintf(stderr, "%s\n", test->get_error_message().c_str());
  }
  test->finalize();
  
  delete test;
  
  if(opt.out) {
    std::ofstream solfile(solfilename.c_str());
    std::stringstream solsstrm;
    solsstrm.precision(std::numeric_limits<VT>::digits10 + 1);
    solsstrm << std::fixed;
    for(GI r = 0; r < matstr->nrows; ++r) {
      spmv_pack_helper::print(solsstrm, x1ptrs, r) << std::endl;
    }
    solfile << solsstrm.str();
    solfile.close();
  }
  
  x1.deallocate();
  x0.deallocate();
  
#else
  std::vector<VT> x1(matstr->nrows), x0(matstr->ncols);
  std::vector<MT> a(matstr->nnz);
  for(GI r = 0; r < matstr->nrows; ++r) {
    x1[r] = 0.0;
  }
  std::mt19937 rng;
  std::uniform_real_distribution<> dis(1.0, 10.0);
  for(GI c = 0; c < matstr->ncols; ++c) {
    x0[c] = dis(rng);
  }
  for(GI i = 0; i < matstr->nnz; ++i) {
    a[i] = 1.0;
  }
  
  spmv_test<GI, MT, VT> * test = nullptr;
  switch(opt.version) {
  case SPMV_TEMPORAL:
    test = new spmv_test_temporal_v1<GI, MT, VT>();
    break;
  case SPMV_QSM0_V1:
    {
      auto obj = new spmv_test_qsm0_v1<GI, LI, DI, MT, VT>();
      obj->set_parameters(
        opt.lookahead, opt.lookback, opt.ntailbits, opt.nbitsperstage
      );
      test = obj;
    }
    break;
  case SPMV_QSM0_V2:
    {
      auto obj = new spmv_test_qsm0_v2<GI, LI, DI, MT, VT>();
      obj->set_parameters(
        opt.lookahead, opt.lookback, opt.ntailbits, opt.nbitsperstage
      );
      test = obj;
    }
    break;
  }
  
  test->init(
    opt.nrepeat, numthreads, &parts[0],
    matstr, &x1[0], &a[0], &x0[0]
  );
  int retval = test->preprocess();
  if(retval == 0) {
#ifdef USING_ITT
    __itt_resume();
#endif
    test->execute();
#ifdef USING_ITT
    __itt_pause();
#endif
    
    test->calculate_metrics();
    
    nflops = test->get_flops_count();
    data_volume = test->get_data_volume();
    mintime = test->get_min_time();
    maxtime = test->get_max_time();
    avgtime = test->get_avg_time();
    flop2byte = test->get_flop2byte();
    bandwidth_mbps = test->get_bandwidth();
    test_complete = true;
    
    test->postprocess();
  } else {
    fprintf(stderr, "%s\n", test->get_error_message().c_str());
  }
  test->finalize();
  
  delete test;
  
  if(opt.out) {
    std::ofstream solfile(solfilename.c_str());
    std::stringstream solsstrm;
    solsstrm.precision(std::numeric_limits<VT>::digits10 + 1);
    solsstrm << std::fixed;
    for(GI r = 0; r < matstr->nrows; ++r) {
      solsstrm << x1[r] << std::endl;
    }
    solfile << solsstrm.str();
    solfile.close();
  }
  
#endif // #ifdef SPMV_PACK_TEST
  
  if(test_complete) {
    printf("Time                : min = %lf s\n", mintime);
    printf("                    : max = %lf s\n", maxtime);
    printf("                    : avg = %lf s\n", avgtime);
    printf("Data volume         : %lf MB\n", data_volume/1024.0/1024.0);
    printf("Effective bandwidth : %lf MB/s\n", bandwidth_mbps);
    printf("Estimated flops     : %lu\n", nflops);
    printf("Estimated flop/byte : %lf\n", flop2byte);
    
    if(opt.batchfilename != nullptr) {
      std::ofstream batchfile;
      std::ifstream tmpbatchfile(opt.batchfilename);
      if(tmpbatchfile.is_open()) {
        tmpbatchfile.close();
        batchfile.open(opt.batchfilename, std::ofstream::out|std::ofstream::app);
      } else {
        batchfile.open(opt.batchfilename, std::ofstream::out);
        batchfile << "casename,version,lookback,lookahead,nbitsperstage,ntailbits,"
                  << "nrepeat,numthreads,OMP_NUM_THREADS,KMP_AFFINITY,"
                  << "mintime,maxtime,avgtime,datavolume,bandwidth,flops,flops/byte"
                  << std::endl;
      }
      
      batchfile << opt.casename << ","
                << version_names[opt.version] << ","
                << opt.lookback << ","
                << opt.lookahead << ","
                << opt.nbitsperstage << ","
                << opt.ntailbits << ","
                << opt.nrepeat << ","
                << numthreads << ","
                << env_omp_num_threads << ","
                << env_kmp_affinity << ","
                << test->get_min_time() << ","
                << test->get_max_time() << ","
                << test->get_avg_time() << ","
                << test->get_data_volume() << ","
                << test->get_bandwidth() << ","
                << test->get_flops_count() << ","
                << test->get_flop2byte() << std::endl;
      
      batchfile.close();
    }
  }
  
  delete matstr;
  
  return 0;
}
