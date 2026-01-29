#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <random>
#include <unistd.h>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/hyb_vector.hpp>

#ifdef USE_OZ
#include <oz_vector.hpp>
#endif

#ifdef USE_NAIVE_HYB
#include "wrapper_naive_hyb.hpp"
#endif

#ifdef USE_ZOMBIT
#include "wrapper_zombit.hpp"
#endif

#ifdef USE_PASTA
#include "wrapper_pasta.hpp"
#endif

#ifdef USE_LA_VECTOR
#include "wrapper_la_vector.hpp"
#endif

#ifndef INDEX_TYPE
#define INDEX_TYPE sdsl::hyb_vector<16>
#endif

#ifdef SEED
std::uint64_t seed = SEED;
#else
std::uint64_t seed = std::random_device{}();
#endif
std::mt19937_64 generator(seed);

// Generate a random number in range [low, high].
std::uint64_t random_number(std::uint64_t low, std::uint64_t high) {
  std::uniform_int_distribution<std::uint64_t> distribution(low, high);
  return distribution(generator);
}

template <typename t_hyb_vec>
void test_random(uint64_t length, int keybit, int density, int n_strings, int n_queries) {
  std::cout << "TESTING: "
    << "length=" << length
    << ", keybit=" << keybit
    << ", density=" << density
    << ", n_strings=" << n_strings
    << ", n_queries=" << n_queries
    << std::endl;
  for (int i_string = 0; i_string < n_strings; ++i_string) {
    // Generate random bitvector.
    sdsl::bit_vector v(length);
    for (uint64_t i = 0; i < length; ++i) {
      v[i] = (random_number(1, density) == 1) ? keybit : 1 - keybit;
    }

    // Initialize hybrid bitvector.
    t_hyb_vec hybrid_v(v);
    typename t_hyb_vec::select_0_type select_0_hybrid_v(&hybrid_v);
    typename t_hyb_vec::select_1_type select_1_hybrid_v(&hybrid_v);

    // Initialize other bitvector for checking.
    sdsl::rrr_vector<15> rrr_v(v);
    sdsl::rrr_vector<15>::rank_0_type rank_0_rrr_v(&rrr_v);
    sdsl::rrr_vector<15>::rank_1_type rank_1_rrr_v(&rrr_v);
    sdsl::rrr_vector<15>::select_0_type select_0_rrr_v(&rrr_v);
    sdsl::rrr_vector<15>::select_1_type select_1_rrr_v(&rrr_v);

    // Compare for select 1 type
    {
      uint64_t one_count = rank_1_rrr_v.rank(length);
      if (one_count > 0) {
        for (int query_i = 0; query_i < n_queries; ++query_i) {
          long i = random_number(1, one_count);
          if (query_i == 0) {
            i =  one_count;
          }
          // for (uint64_t i = 0; i < length; ++i) {
          //   std::cout << "v[" << i << "] = " << v[i] << std::endl;
          // }
          long ret_hybrid = select_1_hybrid_v.select(i);
          long ret_rrr = select_1_rrr_v.select(i);
          if (ret_hybrid != ret_rrr) {
            std::cout << "query_i = " << query_i << std::endl;
            std::cout << "i = " << i << std::endl;
            std::cout << "select = " << "1" << std::endl;
            std::cout << "ret_hybrid = " << ret_hybrid << std::endl;
            std::cout << "ret_rrr = " << ret_rrr << std::endl;
            // sdsl::hyb_vector<16> hybrid_v(v);
            // long ret_hybrid = select_1_hybrid_v.select(i);
            std::exit(EXIT_FAILURE);
          }
        }
      }
    }

    // Compare for select 0 type
    {
      long zero_count = rank_0_rrr_v.rank(length);
      if (zero_count > 0) {
        for (int query_i = 0; query_i < n_queries; ++query_i) {
          long i = random_number(1, zero_count);
          if (query_i == 0) {
            i =  zero_count;
          }
          // for (uint64_t i = 0; i < length; ++i) {
          //   std::cout << "v[" << i << "] = " << v[i] << std::endl;
          // }
          long ret_hybrid = select_0_hybrid_v.select(i);
          long ret_rrr = select_0_rrr_v.select(i);
          if (ret_hybrid != ret_rrr) {
            std::cout << "query_i = " << query_i << std::endl;
            std::cout << "i = " << i << std::endl;
            std::cout << "select = " << "0" << std::endl;
            std::cout << "ret_hybrid = " << ret_hybrid << std::endl;
            std::cout << "ret_rrr = " << ret_rrr << std::endl;
            // sdsl::hyb_vector<16> hybrid_v(v);
            // long ret_hybrid = select_0_hybrid_v.select(i);
            std::exit(EXIT_FAILURE);
          }
        }
      }
    }
  }
}

int main() {
  // Short vectors
  {
    const uint64_t max_length = 1 << 20;
    const uint64_t max_density = 1 << 16;
    const uint64_t n_queries = 500;
    const uint64_t n_strings = 50;

    for (uint64_t length = 1; length <= max_length; length *= 2) {
      for (uint64_t density = 1; density <= max_density; density *= 2) {
        test_random<INDEX_TYPE>(length, 1, density, n_strings, n_queries);
        test_random<INDEX_TYPE>(length, 0, density, n_strings, n_queries);
      }
    }

  }

#ifdef TEST_LARGE_FILES
  // Long vectors
  {
    const uint64_t length = (1ULL << 31) + (1ULL << 10);
    const uint64_t max_density = (1ULL << 20);
    const uint64_t n_queries = 500;
    const uint64_t n_strings = 1;
    for (uint64_t density = 1; density <= max_density; density *= 32) {
      test_random<INDEX_TYPE>(length, 1, density, n_strings, n_queries);
      test_random<INDEX_TYPE>(length, 0, density, n_strings, n_queries);
    }
  }
#endif
}

