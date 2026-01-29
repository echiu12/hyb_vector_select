#include "sdsl/hyb_vector.hpp"
#include "sdsl/io.hpp"
#include <cstdlib>
#include <cstring>
#include <libgen.h>
#include <random>
#include <sys/stat.h>
#include <sys/types.h>

#include <sdsl/bit_vectors.hpp>
#include <sdsl/csa_wt.hpp>

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

#include "benchmark.hpp"

#ifndef VERSION
#define VERSION 0
#endif

#ifndef INDEX_ALIAS
#define INDEX_ALIAS bit_vector
#endif

#ifndef INDEX_TYPE
#define INDEX_TYPE sdsl::bit_vector
#endif

#ifndef BLOCK_SIZE
#define BLOCK_SIZE
#endif

#ifndef NUM_QUERIES
#define NUM_QUERIES 100000
#endif

#ifndef NUM_REPS
#define NUM_REPS 100
#endif

#ifndef CACHE_DIR
#define CACHE_DIR "./cache"
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


// Index for PLCP queries.
template <
  class t_plcp = sdsl::bit_vector,
  typename t_select_1 = typename t_plcp::select_1_type
>
class lcp_index {
public:

  // Aliases.
  using size_type = std::uint64_t;
  using plcp_type = t_plcp;
  using select_support_type = t_select_1;

private:

  // Members.
  plcp_type m_plcp;
  select_support_type m_select_support;

public:

  // Constructor.
  lcp_index(sdsl::cache_config& config) {

    sdsl::bit_vector data;

    // Construct PLCP if not in cache. Otherwise, load PLCP from cache.
    if (!cache_file_exists("plcp", config)) {
      // Construct ISA if not exists.
      if (!cache_file_exists(sdsl::conf::KEY_ISA, config)) {
        sdsl::construct_isa(config);
      }
      sdsl::int_vector_buffer<> isa_buf(cache_file_name(sdsl::conf::KEY_ISA, config));
      // LCP (assume cached from ISA construction)
      sdsl::int_vector<> lcp;
      load_from_file(lcp, cache_file_name(sdsl::conf::KEY_LCP, config));
      // Construct PLCP
      size_type n = lcp.size();
      data = sdsl::bit_vector(2*n, 0);
      size_type data_cnt=0;
      for (size_type i=0, l=0, old_l=1; i < n; ++i) {
          l = lcp[isa_buf[i]];
          data_cnt += l + 1 - old_l;
          data[data_cnt++] = 1;
          old_l = l;
      }
      data.resize(data_cnt);
      sdsl::store_to_cache(data, "plcp", config);
    }
    else {
      sdsl::load_from_cache(data, "plcp", config);
    }

    // Initialize bitvector.
    m_plcp = plcp_type(data);

    // Initialize select support.
    sdsl::util::init_support(m_select_support, &m_plcp);

  }

  // PLCP query.
  std::uint64_t plcp(size_type j) const {
    size_type s = m_select_support.select(j+1);
    return s-(j<<1);
  }

  //! Serialize to a stream.
  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="")const {
      sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
      size_type written_bytes = 0;
      written_bytes += m_plcp.serialize(out, child, "plcp");
      written_bytes += m_select_support.serialize(out, child, "plcp_select");
      sdsl::structure_tree::add_size(child, written_bytes);
      return written_bytes;
  }

  //! Load from a stream.
  void load(std::istream& in) {
      m_plcp.load(in);
      m_select_support.load(in, &m_plcp);
  }

};


int main(int argc, char **argv) {

  // Define parameters.
  auto n_queries = uint64_t{NUM_QUERIES};
  auto n_reps = uint64_t{NUM_REPS};

  // Ensure input file is provided.
  if (argc <= 1) {
    std::cerr << "Error: missing input file" << std::endl;
    exit(EXIT_FAILURE);
  }

  // Get file name.
  auto text_path = std::string(argv[1]);
  char *text_path_raw = strdup(argv[1]);
  auto text_name = std::string(basename(text_path_raw));
  free(text_path_raw);
  std::hash<std::string> hasher;
  auto text_path_hash = std::to_string(hasher(text_path));

  // Get file size.
  auto file_size_bytes = get_file_size_bytes(text_path);

  // Generate queries
  std::vector<std::uint64_t> query_data(n_queries);
  for (std::uint64_t i_query = 0; i_query < n_queries; ++i_query) {
    query_data[i_query] = random_number(0, file_size_bytes - 1);
  }

  // Define cache config.
  mkdir(CACHE_DIR, S_IRWXU);
  sdsl::cache_config config(false, CACHE_DIR, text_path_hash);

  // Populate cache.
  const std::uint8_t width = 8;
  const std::uint8_t num_bytes = 1;
  const char* KEY_TEXT = sdsl::key_text_trait<width>::KEY_TEXT;
  const char* KEY_BWT  = sdsl::key_bwt_trait<width>::KEY_BWT;
  typedef sdsl::int_vector<width> text_type;
  // Load text into cache.
  {
    if (!cache_file_exists(KEY_TEXT, config)) {
      text_type text;
      load_vector_from_file(text, text_path, num_bytes);
      if (contains_no_zero_symbol(text, text_path)) {
        append_zero_symbol(text);
        store_to_cache(text, KEY_TEXT, config);
      }
    }
    register_cache_file(KEY_TEXT, config);
  }
  // Load SA into cache.
  {
    if (!cache_file_exists(sdsl::conf::KEY_SA, config)) {
        sdsl::construct_sa<width>(config);
    }
    register_cache_file(sdsl::conf::KEY_SA, config);
  }
  // Load BWT into cache.
  {
    if (!cache_file_exists(KEY_BWT, config)) {
        sdsl::construct_bwt<width>(config);
    }
    register_cache_file(KEY_BWT, config);
  }
  // Load CSA into cache.
  {
    sdsl::csa_wt<> csa;
    if (!cache_file_exists(std::string(sdsl::conf::KEY_CSA)+"_"+sdsl::util::class_to_hash(csa), config)) {
      construct(csa, text_path, config, num_bytes);
      store_to_cache(csa, std::string(sdsl::conf::KEY_CSA)+"_"+sdsl::util::class_to_hash(csa), config);
    }
    register_cache_file(std::string(sdsl::conf::KEY_CSA)+"_"+sdsl::util::class_to_hash(csa), config);
  }
  // Load LCP into cache.
  {
    if (!cache_file_exists(sdsl::conf::KEY_LCP, config)) {
      construct_lcp_semi_extern_PHI(config);
    }
    register_cache_file(sdsl::conf::KEY_LCP, config);
  }

  // Set namespace.
  using namespace sdsl;

  // Construct index.
  lcp_index<INDEX_TYPE> index(config);

  // Initialize checksum.
  std::uint64_t checksum = 0UL;

  // Warmup
  for (std::uint64_t i_query = 0; i_query < n_queries * 3; ++i_query) {
    std::uint64_t i = query_data[i_query % n_queries];
    std::uint64_t result = index.plcp(i);
    checksum += result;
  }

  // Execute.
  auto start = BENCHMARK_NOW;
  for (std::uint64_t i_query = 0; i_query < n_queries * n_reps; ++i_query) {
    std::uint64_t i = query_data[i_query % n_queries];
    std::uint64_t result = index.plcp(i);
    checksum += result;
  }
  auto end = BENCHMARK_NOW;

  // To prevent compiler optimization.
  if (checksum == 0UL) {
    std::cerr << "Error: Checkum equals zero" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Get results.
  auto version = uint64_t{VERSION};
  auto index_alias = EXPAND_AND_STRINGIFY(INDEX_ALIAS);
  auto index_type = EXPAND_AND_STRINGIFY(INDEX_TYPE);
  auto block_size = EXPAND_AND_STRINGIFY(BLOCK_SIZE);
  sdsl::nullstream ns;
  auto space_bytes = index.serialize(ns);
  auto space_pct = static_cast<double>(space_bytes) / file_size_bytes * 100;
  auto time_total_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count();
  auto time_total_s = time_total_ns / 1e9;
  auto time_us = time_total_ns / (n_reps * n_queries) / 1e3;


  // Store results in map.
  std::map<std::string, std::string> map = {
    { entry("version"),         entry(version) },
    { entry("text_name"),       entry(text_name) },
    { entry("text_path"),       entry(text_path) },
    { entry("index_alias"),     entry(index_alias) },
    { entry("index_type"),      entry(index_type) },
    { entry("block_size"),      entry(block_size) },
    { entry("file_size_bytes"), entry(file_size_bytes) },
    { entry("space_bytes"),     entry(space_bytes) },
    { entry("space_pct"),       entry(space_pct) },
    { entry("time_total_s"),    entry(time_total_s) },
    { entry("time_us"),         entry(time_us) },
  };

  // Print results.
  print_results(map);

}


