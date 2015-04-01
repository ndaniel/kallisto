#ifndef KALLISTO_EMALGORITHM_H
#define KALLISTO_EMALGORITHM_H

#include "common.h"
#include "KmerIndex.h"
#include "weights.h"

#include <algorithm>
#include <iostream>
#include <limits>
#include <vector>


// smallest weight we expect is ~10^-4
// on most machines, TOLERANCE should be 2.22045e-15
//const double TOLERANCE = std::numeric_limits<double>::epsilon() * 10;
//const double TOLERANCE = 1e-100;
const double TOLERANCE = std::numeric_limits<double>::denorm_min();

struct EMAlgorithm {
  // ecmap is the ecmap from KmerIndex
  // counts is vector from collector, with indices corresponding to ec ids
  // target_names is the target_names_ from collector
  // TODO: initialize alpha a bit more intelligently
  EMAlgorithm(const EcMap& ecmap,
              const std::vector<double>& counts,
              const std::vector<std::string>& target_names,
              const std::vector<double>& eff_lens,
              const WeightMap& wm) :
    //idx_(idx),
    num_trans_(target_names.size()),
    ecmap_(ecmap),
    counts_(counts),
    target_names_(target_names),
    eff_lens_(eff_lens),
    weight_map_(wm),
    alpha_(num_trans_, 1.0/num_trans_), // uniform distribution over transcripts
    rho_(num_trans_, 0.0),
    rho_set_(false)
  {
      assert(target_names_.size() == eff_lens.size());
  }

  ~EMAlgorithm() {}

  void run(const std::string& out_dir, size_t n_iter = 500) {
    std::vector<double> next_alpha(alpha_.size(), 0.0);

    assert(weight_map_.size() <= counts_.size());

    double denom;

    std::cout << "[em]\tfishing for the right mixture (. = 50 rounds)" <<
              std::endl;

    std::vector<double> log_lik;
    log_lik.reserve(n_iter);
    std::ofstream out_alpha;
    out_alpha.open(out_dir + "/alpha.txt", std::ios::out);
    if (!out_alpha.is_open()) {
      std::cerr << "Error opening alpha.txt" << std::endl;
      exit(1);
    }

    // output target names to out_alpha
    bool first = true;
    for (auto& t: target_names_) {
      if (!first) {
        out_alpha << "\t" << t;
      } else {
        out_alpha << t;
        first = false;
      }
    }
    out_alpha << std::endl;

    double sum_alpha = 0.0;

    for (auto i = 0; i < n_iter; ++i) {
      if (i % 50 == 0) {
        std::cout << ".";
        std::cout.flush();
        if (i % 500 == 0 && i > 0) {
          std::cout << std::endl;
        }
      }
      sum_alpha = 0.0;
      for (auto a : alpha_) {
        sum_alpha += a;
      }

      double cur_log_lik = 0.0;
      for (auto& ec_kv : ecmap_ ) {
        denom = 0.0;

        // first, compute the denominator: a normalizer
        // iterate over transcripts in EC map
        auto w_search = weight_map_.find(ec_kv.first);

        // everything in ecmap should be in weight_map
        assert( w_search != weight_map_.end() );
        assert( w_search->second.size() == ec_kv.second.size() );

        for (auto t_it = 0; t_it < ec_kv.second.size(); ++t_it) {
          denom += alpha_[ec_kv.second[t_it]] * w_search->second[t_it];
        }

        if (denom < TOLERANCE) {
          continue;
        }

        cur_log_lik += counts_[ec_kv.first] * (log(denom) - log(sum_alpha));

        // compute the update step
        for (auto t_it = 0; t_it < ec_kv.second.size(); ++t_it) {
          next_alpha[ec_kv.second[t_it]] += counts_[ec_kv.first] *
                                            ((w_search->second[t_it] * alpha_[ec_kv.second[t_it]]) / denom);
        }
      }

      // TODO: check for relative difference for convergence in EM
      log_lik.push_back( cur_log_lik );

      // reassign alpha_ to next_alpha
      std::copy(next_alpha.begin(), next_alpha.end(), alpha_.begin());

      write_alpha(out_alpha);

      // clear all next_alpha values 0 for next iteration
      std::fill(next_alpha.begin(), next_alpha.end(), 0.0);
    }

    out_alpha.flush();
    out_alpha.close();
    write_likelihood(log_lik, out_dir + "/ll.txt");

    std::cout << std::endl;
    std::cout.flush();
  }

  void write_likelihood(const std::vector<double>& lik, const std::string& out_fname) {
    std::ofstream out;
    out.open(out_fname, std::ios::out);
    for ( auto l : lik ) {
      out << l << std::endl;
    }

    out.close();
  }
  void compute_rho() {
    if (rho_set_) {
      // rho has already been set, let's clear it
      std::fill(rho_.begin(), rho_.end(), 0.0);
    }

    double total {0.0};
    for (auto i = 0; i < alpha_.size(); ++i) {
      if (eff_lens_[i] < TOLERANCE) {
        std::cerr << "Should actually never really get here... tid: "  << i <<
            std::endl;
        continue;
      }
      rho_[i] = alpha_[i] / eff_lens_[i];
      total += rho_[i];
    }

    for (auto& r : rho_) {
      r /= total;
    }

    rho_set_ = true;
  }

  void write_alpha(std::ofstream& out) {
    bool first = true;
    for (auto a : alpha_) {
      if (!first) {
        out << "\t" << a;
      } else {
        out << a;
        first = false;
      }
    }
    out << std::endl;
  }

  void write(const std::string& out_fname) const {
    std::ofstream out;
    out.open(out_fname, std::ios::out);

    if (!out.is_open()) {
      std::cerr << "Error opening '" << out_fname << "'" <<
          std::endl;
      exit(1);
    }

    out.precision(15);

    out <<
        "target_id" << "\t" <<
        "kallisto_id" << "\t" <<
        "rho" << "\t" <<
        "tpm" << "\t" <<
        "est_counts" <<
        std::endl;

    const double MILLION = 1e6;

    for (auto i = 0; i < rho_.size(); ++i) {
      out <<
          target_names_[i] << "\t" <<
          i << "\t" <<
          rho_[i] << "\t" <<
          rho_[i] * MILLION << "\t" <<
          alpha_[i] <<
          std::endl;
    }

    out.flush();
    out.close();
  }

  int num_trans_;
  const EcMap& ecmap_;
  const std::vector<double>& counts_;
  const std::vector<std::string>& target_names_;
  const std::vector<double>& eff_lens_;
  const WeightMap& weight_map_;
  std::vector<double> alpha_;
  std::vector<double> rho_;
  bool rho_set_;
};

#endif // KALLISTO_EMALGORITHM_H
