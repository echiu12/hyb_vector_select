#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <unistd.h>
#include <cstring>
#include <algorithm>
#include <sys/stat.h>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/csa_wt.hpp>

#define STRINGIFY(x...) #x
#define TYPE_TO_STR(type...) STRINGIFY(type)

#define TEST(t_wt,t_bv...) \
  do { \
    using bv_type = sdsl::t_bv; \
    sdsl::wt_ ## t_wt<bv_type, bv_type::rank_1_type, bv_type::select_1_type, bv_type::select_0_type, sdsl::byte_tree<>> wt; \
    construct(wt, cache_file_name(sdsl::key_trait<width>::KEY_BWT, config), config, 1); \
    test(wt, iteration, fname, length, TYPE_TO_STR(t_wt ## _ ## t_bv), n_queries, query_data, n_reps); \
  } while(0)

#define CACHE_DIR "./cache"

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

// Test lcp queries using given index.
template<typename index_type>
void test(
  index_type &index,
  std::uint64_t iteration,
  std::string fname,
  std::uint64_t length,
  std::string index_name,
  std::uint64_t n_queries,
  std::vector<std::uint64_t> query_data,
  std::uint64_t n_reps
) {
  
  // Get space.
  sdsl::nullstream ns;
  std::uint64_t index_size = sdsl::serialize(index, ns);

  // Log type and space.
  std::cout << std::setw(15) << fname;
  std::cout << std::setw(10) << iteration;
  std::cout << std::setw(40) << index_name;
  std::cout << std::fixed << std::setprecision(4) << std::setw(15) << index_size / 1e6;
  std::cout << std::fixed << std::setprecision(4) << std::setw(15) << (double)index_size / length * 100;
  std::cout << std::flush;
  
  // Initialize checksum.
  std::uint64_t checksum = 0UL;

  // Warmup.
  for (std::uint64_t i_query = 0; i_query < n_queries * 3; ++i_query) {
    uint64_t query_index = 2 * (i_query % n_queries);
    uint64_t c = query_data[query_index];
    uint64_t i = query_data[query_index + 1];
    uint64_t result = index.select(i, c);
    checksum += result;
  }

  // Measure query time.
  auto start = std::chrono::high_resolution_clock::now();
  for (std::uint64_t i_query = 0; i_query < n_queries * n_reps; ++i_query) {
    uint64_t query_index = 2 * (i_query % n_queries);
    uint64_t c = query_data[query_index];
    uint64_t i = query_data[query_index + 1];
    uint64_t result = index.select(i, c);
    checksum += result;
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::uint64_t time_total_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
  auto time_total_s = time_total_ns / 1e9;
  auto time_us = time_total_ns / (n_reps * n_queries) / 1e3;

  // Log times.
  std::cout << std::fixed << std::setprecision(4) << std::setw(15) << time_total_s;
  std::cout << std::fixed << std::setprecision(4) << std::setw(15) << time_us;
  std::cout << std::endl;
  
  // To prevent compiler optimization.
  if (checksum == 0UL)
    std::exit(EXIT_FAILURE);
}

int main(int argc, char **argv) {

  // Define parameters.
  const std::uint64_t n_queries = 100000;
  const std::uint64_t n_reps = 100;
  const std::uint64_t n_iterations = 1;

  // Log header line.
  std::cout << std::setw(15) << "text_name";
  std::cout << std::setw(10) << "iteration";
  std::cout << std::setw(40) << "index_type";
  std::cout << std::setw(15) << "space_mb";
  std::cout << std::setw(15) << "space_pct";
  std::cout << std::setw(15) << "time_total_s";
  std::cout << std::setw(15) << "time_us";
  std::cout << std::endl;

  std::ifstream file_list("files/list");
  std::string fname;
  while (std::getline(file_list, fname)) {
    if (fname.length() > 0) {

      // Get file name.
      std::string fpath = "files/" + fname;

      // Get text length.
      std::ifstream file(fpath);
      file.seekg(0, std::ios::end);
      std::uint64_t length = file.tellg();
      file.close();

      // Initialize cache.
      mkdir(CACHE_DIR, S_IRWXU);
      sdsl::cache_config config(false, CACHE_DIR, fname);

      // Construct reference CSA.
      sdsl::csa_wt<> csa;
      construct(csa, fpath, config, 1);

      // Generate queries.
      std::vector<std::uint64_t> query_data(2 * n_queries);
      for (std::uint64_t i_query = 0; i_query < n_queries; ++i_query) {
        uint64_t comp = random_number(1, csa.sigma - 1);
        uint64_t c = csa.comp2char[comp];
        query_data[2 * i_query] = c;
        uint64_t count = csa.C[comp + 1] - csa.C[comp];
        query_data[2 * i_query + 1] = random_number(1, count);
      }

      // Register cache files.
      const std::uint8_t width = 8;
      const std::uint8_t num_bytes = 1;
      const char* KEY_TEXT = sdsl::key_text_trait<width>::KEY_TEXT;
      const char* KEY_BWT  = sdsl::key_bwt_trait<width>::KEY_BWT;
      typedef sdsl::int_vector<width> text_type;
      {
        if (!cache_file_exists(KEY_TEXT, config)) {
          text_type text;
          load_vector_from_file(text, fpath, num_bytes);
          if (contains_no_zero_symbol(text, fpath)) {
            append_zero_symbol(text);
            store_to_cache(text, KEY_TEXT, config);
          }
        }
        register_cache_file(KEY_TEXT, config);
      }
      {
        if (!cache_file_exists(sdsl::conf::KEY_SA, config)) {
            sdsl::construct_sa<width>(config);
        }
        register_cache_file(sdsl::conf::KEY_SA, config);
      }
      {
        if (!cache_file_exists(KEY_BWT, config)) {
            sdsl::construct_bwt<width>(config);
        }
        register_cache_file(KEY_BWT, config);
      }

      for (std::uint64_t iteration = 0; iteration < n_iterations; ++iteration) {

        TEST(blcd, rrr_vector<15>);
        TEST(blcd, rrr_vector<31>);
        TEST(blcd, rrr_vector<63>);
        TEST(blcd, rrr_vector<127>);
        TEST(blcd, hyb_vector<8>);
        TEST(blcd, hyb_vector<16>);
        TEST(blcd, hyb_vector<32>);
        TEST(blcd, hyb_vector<64>);
        TEST(blcd, sd_vector<>);
        TEST(blcd, bit_vector);

        TEST(huff, rrr_vector<15>);
        TEST(huff, rrr_vector<31>);
        TEST(huff, rrr_vector<63>);
        TEST(huff, rrr_vector<127>);
        TEST(huff, hyb_vector<8>);
        TEST(huff, hyb_vector<16>);
        TEST(huff, hyb_vector<32>);
        TEST(huff, hyb_vector<64>);
        TEST(huff, sd_vector<>);
        TEST(huff, bit_vector);

      }

    }

  }

  file_list.close(); // close file list
}

