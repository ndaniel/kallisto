#ifndef KALLISTO_MINCOLLECTOR_H
#define KALLISTO_MINCOLLECTOR_H

#include "common.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <algorithm>


template <typename Index>
struct MinCollector {

  MinCollector(Index& ind, const ProgramOptions& opt) : index(ind), counts(index.ecmap.size(), 0) {}




  void collect(std::vector<int>& v, kseq_t *s1, kseq_t *s2, bool paired) {
    if (v.empty()) {
      return;
    }
    sort(v.begin(), v.end()); // sort by increasing order

    int count = 1; // how many k-mer support the ec
    std::vector<int> u = index.ecmap[v[0]];

    for (int i = 1; i < v.size(); i++) {
      if (v[i] != v[i-1]) {
        u = index.intersect(v[i],u);
        if (u.empty()) {
          break;
        }
      }
      count++; // increase the count
    }
    // if u is empty do nothing
    if (u.empty()) {
      int m1=1,m2=1;
      int c = 1;
      for (int i = 1; i < v.size(); i++) {
        if (v[i] != v[i-1]) {
          c = 1;
        } else {
          c++;
          if (c > m2) {
            if (c < m1) {
              m2 = c;
            } else {
              m1 = c;
            }
          }
        }
      }

      // only output if two classes have more then 1 kmer supporting
      if (m1 > 1 && m2 > 1) {
        c = 1;
        for (int i = 1; i < v.size(); i++) {
          if (v[i] != v[i-1]) {
            std::cout << "(" << v[i-1] << "," << c << "),";
            c = 1;
          } else {
            c++;
          }
        }
        
        std::cout << "(" << v[v.size()-1] <<"," << c << ")\n";
        std::cout << "@" << s1->name.s << "\n" << s1->seq.s << "\n+\n" << s1->qual.s << "\n";
        if (paired) {
          std::cout << "@" << s2->name.s << "\n" << s2->seq.s << "\n+\n" << s2->qual.s << "\n";
        }
      }
    

      return;
    }

    auto search = index.ecmapinv.find(u);
    if (search != index.ecmapinv.end()) {
      // ec class already exists, update count
      ++counts[search->second];
    } else {
      // new ec class, update the index and count
      auto necs = counts.size();
      index.ecmap.insert({necs,u});
      index.ecmapinv.insert({u,necs});
      counts.push_back(1);
    }
  }

  void write(std::ostream& o) {
    for (int id = 0; id < counts.size(); id++) {
      o << id << "\t" << counts[id] << "\n";
    }
  }

  void loadCounts(ProgramOptions& opt) {
    int num_ecs = counts.size();
    counts.clear();
    std::ifstream in((opt.output + "/counts.txt"));
    int i = 0;
    if (in.is_open()) {
      std::string line;
      while (getline(in, line)) {
        std::stringstream ss(line);
        int j,c;
        ss >> j;
        ss >> c;
        if (j != i) {
          std::cerr << "Error: equivalence class does not match index. Found "
                    << j << ", expected " << i << std::endl;
          exit(1);
        }
        counts.push_back(c);
        i++;
      }

      if (i != num_ecs) {
        std::cerr << "Error: number of equivalence classes does not match index. Found "
                  << i << ", expected " << num_ecs << std::endl;
        exit(1);
      }
    } else {
      std::cerr << "Error: Could not open file " << opt.output << "/counts.txt" << std::endl;
      exit(1);

    }
  }

  Index& index;
  std::vector<int> counts;

};

#endif // KALLISTO_MINCOLLECTOR_H
