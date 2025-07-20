#include <chrono>
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

#define TEST(bv...) \
  do { \
    lcp_index<sdsl::csa_wt<>, sdsl::bv> lcp(config); \
    test(lcp, iteration, fname, length, TYPE_TO_STR(bv), n_queries, query_data, n_reps); \
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

// Index for PLCP queries.
template <class t_csa = sdsl::csa_wt<>, 
          class t_plcp = sdsl::bit_vector, 
          class t_select_support = typename t_plcp::select_1_type>
class lcp_index {
public:
  using size_type = std::uint64_t;
private:
#ifndef USE_PLCP
  t_csa m_csa;
#endif
  t_plcp m_plcp;
  t_select_support m_select_support;
public:
  lcp_index(sdsl::cache_config& config) {
    sdsl::bit_vector data;
#ifndef USE_PLCP
    sdsl::load_from_cache(m_csa, std::string(sdsl::conf::KEY_CSA)+"_"+sdsl::util::class_to_hash(m_csa), config);
#endif
    if (!cache_file_exists("plcp", config)) {
      // ISA
      if (!cache_file_exists(sdsl::conf::KEY_ISA, config)) {
        construct_isa(config);
      }
      sdsl::int_vector_buffer<> isa_buf(cache_file_name(sdsl::conf::KEY_ISA, config));
      // LCP (assume cached)
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
    m_plcp = t_plcp(data);
    sdsl::util::init_support(m_select_support, &m_plcp);
  }

#ifndef USE_PLCP
    std::uint64_t lcp(size_type i) const {
        size_type j = m_csa[i];
        size_type s = m_select_support.select(j+1);
        return s-(j<<1);
    }
#endif

  std::uint64_t plcp(size_type j) const {
    size_type s = m_select_support.select(j+1);
    return s-(j<<1);
  }

  //! Serialize to a stream.
  size_type
  serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr,
            std::string name="")const {
      sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name,
                                    sdsl::util::class_name(*this));
      size_type written_bytes = 0;
      written_bytes += m_plcp.serialize(out, child, "plcp");
      written_bytes += m_select_support.serialize(out, child,
                        "plcp_select");
      sdsl::structure_tree::add_size(child, written_bytes);
      return written_bytes;
  }

  //! Load from a stream.
  void load(std::istream& in, const t_csa* csa) {
#ifndef USE_PLCP
      m_csa = csa;
#endif
      m_plcp.load(in);
      m_select_support.load(in, &m_plcp);
  }

};

// Test LCP/PLCP queries.
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
  
  // Warmup
  for (std::uint64_t i_query = 0; i_query < n_queries * 3; ++i_query) {
    std::uint64_t i = query_data[i_query % n_queries];
#ifdef USE_PLCP
    std::uint64_t result = index.plcp(i);
#else 
    std::uint64_t result = index.lcp(i);
#endif // USE_PLCP
    checksum += result;
  }

  // Measure query time.
  auto start = std::chrono::high_resolution_clock::now();
  for (std::uint64_t i_query = 0; i_query < n_queries * n_reps; ++i_query) {
    std::uint64_t i = query_data[i_query % n_queries];
#ifdef USE_PLCP
    std::uint64_t result = index.plcp(i);
#else 
    std::uint64_t result = index.lcp(i);
#endif // USE_PLCP
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

  // Define parameters.s
#ifdef USE_PLCP
  const std::uint64_t n_queries = 100000;
  const std::uint64_t n_reps = 100;
  const std::uint64_t n_iterations = 1;
#else
  const std::uint64_t n_queries = 100000;
  const std::uint64_t n_reps = 100;
  const std::uint64_t n_iterations = 1;
#endif // USE_PLCP

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

      // Get file path.
      std::string fpath = "files/" + fname;

      // Get text length.
      std::ifstream file(fpath);
      file.seekg(0, std::ios::end);
      std::uint64_t length = file.tellg();
      file.close();

      // Generate queries
      std::vector<std::uint64_t> query_data(n_queries);
      for (std::uint64_t i_query = 0; i_query < n_queries; ++i_query) {
        query_data[i_query] = random_number(0, length - 1);
      }

      // Initialize cache
      mkdir(CACHE_DIR, S_IRWXU);
      sdsl::cache_config config(false, CACHE_DIR, fname);
      
      // Register cache files
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
      {
        sdsl::csa_wt<> csa;
        if (!cache_file_exists(std::string(sdsl::conf::KEY_CSA)+"_"+sdsl::util::class_to_hash(csa), config)) {
          construct(csa, fpath, config, num_bytes);
          store_to_cache(csa, std::string(sdsl::conf::KEY_CSA)+"_"+sdsl::util::class_to_hash(csa), config);
        }
        register_cache_file(std::string(sdsl::conf::KEY_CSA)+"_"+sdsl::util::class_to_hash(csa), config);
      }
      {
        if (!cache_file_exists(sdsl::conf::KEY_LCP, config)) {
          construct_lcp_semi_extern_PHI(config);
        }
        register_cache_file(sdsl::conf::KEY_LCP, config);
      }

      for (std::uint64_t iteration = 0; iteration < n_iterations; ++iteration)
      {
        TEST(rrr_vector<15>);
        TEST(rrr_vector<31>);
        TEST(rrr_vector<63>);
        TEST(rrr_vector<127>);
        TEST(hyb_vector<8>);
        TEST(hyb_vector<16>);
        TEST(hyb_vector<32>);
        TEST(hyb_vector<64>);
        TEST(sd_vector<>);
        TEST(bit_vector);
      }

    }

  }

  file_list.close(); // close file list
}

