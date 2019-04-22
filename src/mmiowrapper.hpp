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

#include <cstdio>
#include "mmio.h"

#include <map>

enum mm_object {
  MMO_MATRIX = 0,
  MMO_OTHER
};

enum mm_format {
  MMF_COORDINATE = 0,
  MMF_ARRAY
};

enum mm_type {
  MMT_REAL = 0,
  MMT_COMPLEX,
  MMT_INTEGER,
  MMT_PATTERN
};

enum mm_symmetry {
  MMS_GENERAL = 0,
  MMS_SYMMETRIC,
  MMS_SKEWSYMMETRIC,
  MMS_HERMITIAN
};

char const * mm_symmetry_names[] = {
  "general",
  "symmetric",
  "skew symmetric",
  "hermitian"
};

struct mm_metadata {
  mm_object object;
  mm_format format;
  mm_type field;
  mm_symmetry symmetry;
  int nrows;
  int ncols;
  int nnz;
};

int read_mm_matrix_metadata(FILE * f, mm_metadata * md) {
  MM_typecode matcode;
  int nrows, ncols, nnz;
  
  if(mm_read_banner(f, &matcode) != 0) {
    return 1; // could not process MatrixMarket banner
  }
  
  if(mm_read_mtx_crd_size(f, &nrows, &ncols, &nnz) != 0) {
    return 2; // could not read matrix dimensions
  }
  
  if(mm_is_matrix(matcode)) {
    md->object = MMO_MATRIX;
  } else {
    md->object = MMO_OTHER;
  }
  
  if(mm_is_coordinate(matcode)) {
    md->format = MMF_COORDINATE;
  } else {
    md->format = MMF_ARRAY;
  }
  
  if(mm_is_real(matcode)) {
    md->field = MMT_REAL;
  } else if(mm_is_complex(matcode)) {
    md->field = MMT_COMPLEX;
  } else if(mm_is_integer(matcode)) {
    md->field = MMT_INTEGER;
  } else if(mm_is_pattern(matcode)) {
    md->field = MMT_PATTERN;
  }
  
  if(mm_is_general(matcode)) {
    md->symmetry = MMS_GENERAL;
  } else if(mm_is_symmetric(matcode)) {
    md->symmetry = MMS_SYMMETRIC;
  } else if(mm_is_skew(matcode)) {
    md->symmetry = MMS_SKEWSYMMETRIC;
  } else if(mm_is_hermitian(matcode)) {
    md->symmetry = MMS_HERMITIAN;
  }
  
  md->nrows = nrows;
  md->ncols = ncols;
  md->nnz = nnz;
  
  return 0;
}

template<typename GI>
int read_mm_matrix_structure(
  FILE * f, mm_metadata const & md, std::pair<GI, GI> * s
) {
  if(md.object != MMO_MATRIX && md.format != MMF_COORDINATE) {
    return 1;
  }
  
  int r, c;
  switch(md.field) {
    case MMT_REAL: {
      double val;
      for(int i = 0; i < md.nnz; ++i) {
        fscanf(f, "%d %d %lf\n", &r, &c, &val); --r; --c;
        s[i] = std::pair<GI, GI>(r, c);
      }
    } break;
    case MMT_COMPLEX: {
      double rl, cm;
      for(int i = 0; i < md.nnz; ++i) {
        fscanf(f, "%d %d %lf %lf\n", &r, &c, &rl, &cm); --r; --c;
        s[i] = std::pair<GI, GI>(r, c);
      }
    } break;
    case MMT_INTEGER: {
      int val;
      for(int i = 0; i < md.nnz; ++i) {
        fscanf(f, "%d %d %d\n", &r, &c, &val); --r; --c;
        s[i] = std::pair<GI, GI>(r, c);
      }
    } break;
    case MMT_PATTERN: {
      for(int i = 0; i < md.nnz; ++i) {
        fscanf(f, "%d %d\n", &r, &c); --r; --c;
        s[i] = std::pair<GI, GI>(r, c);
      }
    } break;
  }
  
  return 0;
}
