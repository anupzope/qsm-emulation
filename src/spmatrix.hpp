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

#include <map>
#include <queue>
#include <limits>
#include <algorithm>

template<typename GI>
struct matrixentries {
  virtual ~matrixentries() = 0;
  virtual void allocate(GI const nnz) = 0;
  virtual void deallocate() = 0;
};

template<typename GI, typename VT>
struct matrixentries_t : matrixentries<GI> {
  VT * values;
  
  matrixentries_t() : values(nullptr) {
  }
  
  ~matrixentries_t() {
    if(values) delete[] values;
    values = nullptr;
  }
  
  void allocate(int const nnz) {
    values = new VT[nnz];
  }
  
  void deallocate() {
    if(values) delete[] values;
    values = nullptr;
  }
};

template<typename GI>
struct spmatrixstructure_crs {
  GI nnz;
  GI nrows;
  GI ncols;
  GI * rowoffs;
  GI * colidxs;
  
  spmatrixstructure_crs() : nnz(0), nrows(0), ncols(0), rowoffs(nullptr),
  colidxs(nullptr) {
  }
  
  ~spmatrixstructure_crs() {
    deallocate();
  }
  
  void deallocate() {
    if(rowoffs) delete[] rowoffs;
    if(colidxs) delete[] colidxs;
    nnz = 0;
    nrows = 0;
    ncols = 0;
    rowoffs = nullptr;
    colidxs = nullptr;
  }
  
  void generate(
    GI const nnz, GI const nrows, GI const ncols,
    std::pair<GI, GI> const * coordinates
  ) {
    deallocate();
    
    this->nnz = nnz;
    this->nrows = nrows;
    this->ncols = ncols;
    
    this->rowoffs = new GI[nrows+1];
    for(GI r = 0; r < nrows+1; ++r) {
      this->rowoffs[r] = 0;
    }
    
    for(GI i = 0; i < nnz; ++i) {
      this->rowoffs[coordinates[i].first+1]++;
    }
    
    for(GI r = 0; r < nrows; ++r) {
      this->rowoffs[r+1] += this->rowoffs[r];
    }
    
    GI * count = new GI[nrows];
    for(GI r = 0; r < nrows; ++r) {
      count[r] = 0;
    }
    
    this->colidxs = new GI[this->rowoffs[nrows]];
    for(GI i = 0; i < nnz; ++i) {
      GI r = coordinates[i].first;
      this->colidxs[this->rowoffs[r]+count[r]++] = coordinates[i].second;
    }
    
    delete[] count;
  }
};

template<typename GI>
void rcmorder(
  GI const nnodes, GI const * nbroffs, GI const * nbridxs, GI * o2n//, GI * n2o
) {
  struct util {
    static void add_unused_neighbors_to_queue(
      GI const n, GI const * nbroffs, GI const * nbridxs, GI const * o2n,
      std::queue<GI> & q
    ) {
      struct cmp {
        bool operator()(
          std::pair<GI, GI> const & l,
          std::pair<GI, GI> const & r
        ) {
          return l.first < r.first;
        }
      };
      
      std::vector<std::pair<GI, GI>> nbrs;
      for(GI i = nbroffs[n]; i < nbroffs[n+1]; ++i) {
        if(o2n[nbridxs[i]] == -1) {
          nbrs.push_back(
            std::pair<GI, GI>(
              nbroffs[nbridxs[i]+1] - nbroffs[nbridxs[i]], nbridxs[i]
            )
          );
        }
      }
      std::sort(nbrs.begin(), nbrs.end(), cmp());
      for(auto const & a : nbrs) {
        q.push(a.second);
      }
    }
  };
  
  std::vector<std::pair<GI, GI>> sorteddeg(nnodes);
  GI sdegidx = 0;
  #pragma omp parallel for
  for(GI n = 0; n < nnodes; ++n) {
    sorteddeg[n] = std::pair<GI, GI>(nbroffs[n+1]-nbroffs[n], n);
  }
  std::sort(sorteddeg.begin(), sorteddeg.end());
  
  std::queue<GI> q;
  GI idx = 0;
  
  #pragma omp parallel for
  for(GI n = 0; n < nnodes; ++n) {
    o2n[n] = -1;
  }
  
  constexpr int nslabs = 100;
  std::vector<GI> slabs(nslabs+1);
  slabs[0] = 0;
  for(int i = 1; i < nslabs+1; ++i) {
    slabs[i] = slabs[i-1] + (nnodes/nslabs) + (i <= nnodes%nslabs ? 1 : 0);
  }
  int slabidx = 0;
  
  bool success = true;
  printf("  Performing Cuthill-McKee...\n");
  while(idx < nnodes) {
    // Select unused lowest degree node as starting node
    GI sn;
    while(sdegidx < nnodes) {
      if(o2n[sorteddeg[sdegidx].second] == -1) {
        sn = sorteddeg[sdegidx].second;
        break;
      }
      ++sdegidx;
    }
    
    o2n[sn] = idx++;
    util::add_unused_neighbors_to_queue(sn, nbroffs, nbridxs, &o2n[0], q);
    while(slabidx <= nslabs && idx >= slabs[slabidx]) {
      printf("\r  %d%%", slabidx);
      fflush(stdout);
      ++slabidx;
    }
    
    while(!q.empty()) {
      GI n = q.front();
      q.pop();
      if(o2n[n] == -1) {
        o2n[n] = idx++;
        util::add_unused_neighbors_to_queue(n, nbroffs, nbridxs, &o2n[0], q);
        while(slabidx <= nslabs && idx >= slabs[slabidx]) {
          printf("\r  %d%%", slabidx);
          fflush(stdout);
          ++slabidx;
        }
      }
    }
  }
  printf("\r  %d%%\n", slabidx-1);
  fflush(stdout);
  
  printf("  Reversing CM sequence...\n");
  #pragma omp parallel for
  for(GI n = 0; n < nnodes; ++n) {
    GI nn = o2n[n];
    nn = nnodes-1-nn;
    o2n[n] = nn;
    //n2o[nn] = n;
  }
}

template<typename GI>
void rcmorder(
  spmatrixstructure_crs<GI> * matstr,
  GI * o2n,
  GI * n2o
) {
  std::queue<GI> q;
  GI idx = 0;
  
  for(GI r = 0; r < matstr->nrows; ++r) {
    o2n[r] = -1;
    n2o[r] = -1;
  }
  
  struct util {
    static GI get_lowest_degree_unused_row(
      GI const nrows, GI const * rowoffs, GI const * o2n
    ) {
      GI mindeg = std::numeric_limits<GI>::max();
      GI mindegrow = -1;
      
      for(GI r = 0; r < nrows; ++r) {
        if(o2n[r] == -1) {
          GI deg = rowoffs[r+1]-rowoffs[r];
          if(mindeg > deg) {
            mindeg = deg;
            mindegrow = r;
          }
        }
      }
      
      return mindegrow;
    }
    
    static void add_unused_neighbors_to_queue(
      GI const row, GI const * rowoffs, GI const * colidxs, GI const * o2n,
      std::queue<GI> & q
    ) {
      std::vector<std::pair<GI, GI>> nbrs;
      for(GI i = rowoffs[row]; i < rowoffs[row+1]; ++i) {
        if(o2n[colidxs[i]] == -1) {
          nbrs.push_back(
            std::pair<GI, GI>(
              rowoffs[colidxs[i]+1] - rowoffs[colidxs[i]], colidxs[i]
            )
          );
        }
      }
      std::sort(nbrs.begin(), nbrs.end());
      for(auto a : nbrs) {
        q.push(a.second);
      }
    }
  } ut;
  
  while(idx < matstr->nrows) {
    GI n = ut.get_lowest_degree_unused_row(matstr->nrows, matstr->rowoffs, &o2n[0]);
    if(n != -1) {
      GI newidx = (matstr->nrows-1)-idx++;
      o2n[n] = newidx;
      n2o[newidx] = n;
      ut.add_unused_neighbors_to_queue(n, matstr->rowoffs, matstr->colidxs, &o2n[0], q);
    }
    
    while(!q.empty()) {
      GI row = q.front();
      q.pop();
      if(o2n[row] == -1) {
        GI newidx = (matstr->nrows-1)-idx++;
        o2n[row] = newidx;
        n2o[newidx] = row;
        ut.add_unused_neighbors_to_queue(row, matstr->rowoffs, matstr->colidxs, &o2n[0], q);
      }
    }
  }
}

template<typename GI>
void gorder(
  GI const nnodes, GI const * inoffs, GI const * inidxs,
  GI const * outoffs, GI const * outidxs, GI * o2n, GI * n2o
) {
  /*struct key {
    GI score;
    GI node;
  };
  
  struct key_comparator {
    bool operator()(key const & l, key const & r) {
      return l.score < r.score;
    }
  };
  
  for(GI n = 0; n < nnodes; ++n) {
    o2n[n] = -1;
    n2o[n] = -1;
  }
  
  // find node with largest in-degree
  GI maxindeg = 0, startnode;
  for(GI n = 0; n < nnodes; ++n) {
    GI indeg = inoffs[n+1]-inoffs[n];
    if(maxindeg < indeg) {
      maxindeg = indeg;
      startnode = n;
    }
  }
  
  std::vector<key> pq;
  for(GI n = 0; n < nnodes; ++n) {
    if(n != startnode) {
      q.push(key(0, n));
    }
  }
  
  std::make_heap(pq.begin(), pq.end());
  GI idx = 0;
  o2n[startnode] = idx;
  n2o[idx] = startnode;
  ++idx;
  while(idx < nnodes) {
    GI ve = n2o[idx-1];
    for(GI i = outoffs[ve]; i < outoffs[ve+1]; ++i) {
      if(o2n[outidxs[i]] == -1) {
        
      }
    }
    for(GI 
  }*/
}

template<typename GI>
void reorder(GI const * o2n, spmatrixstructure_crs<GI> * matstr) {
  GI * rowoffs = new GI[matstr->nrows];
  rowoffs[0] = 0;
  for(GI r = 0; r < matstr->nrows; ++r) {
    GI nr = o2n[r];
    rowoffs[nr+1] = matstr->rowoffs[r+1] - matstr->rowoffs[r];
  }
  
  GI * count = new GI[matstr->nrows];
  for(GI r = 0; r < matstr->nrows; ++r) {
    rowoffs[r+1] += rowoffs[r];
    count[r] = 0;
  }
  
  GI * colidxs = new GI[rowoffs[matstr->nrows]];
  for(GI r = 0; r < matstr->nrows; ++r) {
    GI nr = o2n[r];
    for(GI i = matstr->rowoffs[r]; i < matstr->rowoffs[r+1]; ++i) {
      colidxs[rowoffs[nr]+count[nr]++] = o2n[matstr->colidxs[i]];
    }
  }
  
  delete[] matstr->rowoffs;
  matstr->rowoffs = rowoffs;
  
  delete[] matstr->colidxs;
  matstr->colidxs = colidxs;
  
  delete[] count;
}

template<typename GI>
void calculate_bandwidth(
  GI const nrows, GI const nnz, std::pair<GI, GI> const * coord,
  GI * lobw, GI * upbw, double * loavgbw, double * upavgbw,
  GI * maxdegree, GI * mindegree, double * avgdegree
) {
  std::vector<std::pair<GI, GI>> loupbw(nrows);
  std::vector<GI> degree(nrows);
  
  #pragma omp parallel for
  for(GI r = 0; r < nrows; ++r) {
    loupbw[r] = std::pair<GI, GI>(0, 0);
    degree[r] = 0;
  }
  
  for(GI i = 0; i < nnz; ++i) {
    GI r = coord[i].first;
    GI c = coord[i].second;
    if(r > c) {
      loupbw[r].first = std::max(r-c, loupbw[r].first);
    } else if(r < c) {
      loupbw[r].second = std::max(c-r, loupbw[r].second);
    }
    degree[r]++;
  }
  
  *maxdegree = 0;
  *mindegree = std::numeric_limits<GI>::max();
  *avgdegree = 0.0;
  for(GI i = 0; i < nrows; ++i) {
    *maxdegree = std::max(degree[i], *maxdegree);
    *mindegree = std::min(degree[i], *mindegree);
    *avgdegree += (double)degree[i];
  }
  *avgdegree /= (double)nrows;
  
  *lobw = 0;
  *upbw = 0;
  *loavgbw = 0.0;
  *upavgbw = 0.0;
  for(GI r = 0; r < nrows; ++r) {
    *lobw = std::max(*lobw, loupbw[r].first);
    *upbw = std::max(*upbw, loupbw[r].second);
    *loavgbw += (double)(loupbw[r].first);
    *upavgbw += (double)(loupbw[r].second);
  }
  *loavgbw /= (double)nrows;
  *upavgbw /= (double)nrows;
}
