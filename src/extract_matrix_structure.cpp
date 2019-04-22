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

#include "mmiowrapper.hpp"
#include "spmatrix.hpp"

#include <set>
#include <cstring>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <fstream>

#include <argp.h>

template<typename GI>
int extract_coordinates_from_matrix_market(
  char const * filename,
  std::vector<std::pair<GI, GI>> * coord,
  GI * nrows, GI * ncols, GI * nnz,
  mm_symmetry * symmetry
) {
  coord->clear();
  
  mm_metadata md;
  FILE *f;
  int retval;
  
  printf("Opening matrix file...\n");
  if((f = fopen(filename, "r")) == NULL) {
    fprintf(stderr, "could not open Matrix Market file: %s\n", filename);
    return 1;
  }
  
  printf("Reading matrix metadata...\n");
  if((retval = read_mm_matrix_metadata(f, &md)) != 0) {
    fclose(f);
    switch(retval) {
    case 1:
      fprintf(stderr, "could not read Matrix Market banner\n");
      break;
    case 2:
      fprintf(stderr, "could not read matrix dimensions\n");
      break;
    }
    return 1;
  }
  
  coord->resize(md.nnz);
  
  printf("Reading matrix coordinates...\n");
  if(read_mm_matrix_structure(f, md, &(*coord)[0]) != 0) {
    fprintf(stderr, "could not read matrix coordinates\n");
    fclose(f);
    return 1;
  }
  
  fclose(f);
  
  printf("Generalizing matrix structure...\n");
  switch(md.symmetry) {
  case MMS_SYMMETRIC:
  case MMS_SKEWSYMMETRIC:
  case MMS_HERMITIAN:
    for(GI i = 0; i < md.nnz; ++i) {
      if((*coord)[i].first != (*coord)[i].second) {
        coord->push_back(
          std::pair<GI, GI>((*coord)[i].second, (*coord)[i].first)
        );
      }
    }
    break;
  }
  
  *nrows = md.nrows;
  *ncols = md.ncols;
  *nnz = coord->size();
  *symmetry = md.symmetry;
  
  return 0;
}

template<typename GI>
void write_crs(
  char const * matrix_filename,
  GI const nrows, GI const ncols, GI const nnz,
  std::pair<GI, GI> const * coord
) {
  spmatrixstructure_crs<GI> * matstr = new spmatrixstructure_crs<GI>();
  matstr->generate(nnz, nrows, ncols, &coord[0]);
  
  printf("Writing matrix structure in CRS format...\n");
  //format:
  //nrows ncols nnz
  //n1 colidx1 colidx2...
  //n2 colidx2 colidx2...
  //...
  std::ofstream outfile(matrix_filename);
  std::stringstream outsstrm;
  outfile << matstr->nrows << " " << matstr->ncols << " " << matstr->nnz << std::endl;
  for(GI r = 0; r < matstr->nrows; ++r) {
    GI nconns = matstr->rowoffs[r+1]-matstr->rowoffs[r];
    outsstrm << nconns;
    for(GI i = 0; i < nconns; ++i) {
      outsstrm << " " << matstr->colidxs[matstr->rowoffs[r]+i];
    }
    outsstrm << std::endl;
  }
  outfile << outsstrm.str();
  outfile.close();
  
  delete matstr;
}

struct program_options {
  char const * casename;
  char const * batchfilename;
  
  program_options() {
    casename = nullptr;
    batchfilename = nullptr;
  }
};

static int parse_options(int key, char * arg, argp_state * state) {
  program_options * opt = (program_options*)state->input;
  switch(key) {
  case 403:
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

int main(int argc, char *argv[]) {
  using GI = int;
  
  program_options opt;
  
  struct argp_option options[] = {
    {
      "batch", 403, "STRING", 0,
      "file for appending matrix information"
    },
    {0}
  };
  
  argp a = {options, parse_options, "CASE-NAME"};
  
  int parserv = 0;
  if((parserv = argp_parse(&a, argc, argv, 0, 0, &opt)))
    return parserv;
  
  printf("case name  : %s\n", opt.casename);
  
  std::vector<std::pair<GI, GI>> coord;
  GI nrows, ncols, nnz;
  mm_symmetry symmetry;
  
  std::string infilename(opt.casename);
  infilename += ".mtx";
  
  if(extract_coordinates_from_matrix_market(infilename.c_str(), &coord, &nrows,
    &ncols, &nnz, &symmetry) != 0) {
        fprintf(stderr, "could not extract coordinates from input file\n");
        return 1;
  }
  
  printf("nrows      : %d\n", nrows);
  printf("ncols      : %d\n", ncols);
  printf("nnz        : %d\n", nnz);
  printf("symmetry   : %s\n", mm_symmetry_names[symmetry]);
  printf("\n");
  
  std::sort(coord.begin(), coord.end());
  GI nndnz = 0; // number of non-diagonal non-zeros
  for(GI i = 0; i < nnz; ++i) {
    if(coord[i].first != coord[i].second)
      ++nndnz;
  }
  
  printf("Calculating matrix properties...\n");
  GI lobw, upbw, mindegree, maxdegree;
  double loavgbw, upavgbw, avgdegree;
  calculate_bandwidth(
    nrows, nnz, &coord[0], &lobw, &upbw, &loavgbw, &upavgbw,
    &maxdegree, &mindegree, &avgdegree
  );
  printf("  lower bandwidth : %d, avg = %lf\n", lobw, loavgbw);
  printf("  upper bandwidth : %d, avg = %lf\n", upbw, upavgbw);
  printf("  degree          : min = %d, max = %d, avg = %lf\n",
    mindegree, maxdegree, avgdegree);
  printf("\n");
  
  if(opt.batchfilename != nullptr) {
    std::ofstream batchfile;
    std::ifstream tmpbatchfile(opt.batchfilename);
    if(tmpbatchfile) {
      tmpbatchfile.close();
      batchfile.open(opt.batchfilename, std::ofstream::out|std::ofstream::app);
    } else {
      batchfile.open(opt.batchfilename, std::ofstream::out);
      batchfile << "casename,nrows,ncols,nnz,"
        "lobw,upbw,avg lobw,avg upbw,max degree,min degree,avg degree"
        << std::endl;
    }
    batchfile << opt.casename << "," << nrows << "," << ncols << "," << nnz
              << "," << lobw << "," << upbw << "," << loavgbw << "," << upavgbw
              << "," << maxdegree << "," << mindegree << "," << avgdegree
              << std::endl;
    batchfile.close();
  }
  
  std::string outfilename(opt.casename);
  outfilename += ".crs";
  write_crs(outfilename.c_str(), nrows, ncols, nnz, &coord[0]);
  
  return 0;
}
