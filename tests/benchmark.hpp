#pragma once

#include <cassert>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <map>
#include <vector>

#define STRINGIFY(...) #__VA_ARGS__
#define EXPAND_AND_STRINGIFY(...) STRINGIFY(__VA_ARGS__)

#define REMOVE_FIRST(X,...) __VA_ARGS__
#define EXPAND_AND_REMOVE_FIRST(...) REMOVE_FIRST(__VA_ARGS__)

#define BENCHMARK_NOW std::chrono::high_resolution_clock::now()

#ifndef BENCHMARK_FIELDS
#define BENCHMARK_FIELDS \
  X(text_name) \
  X(index_alias) \
  X(block_size) \
  X(time_us) \
  X(space_pct) \
  X(time_total_s) \
  X(space_bytes) \
  X(file_size_bytes) \
  X(index_type)
#endif // BENCHMARK_FIELDS

uint64_t get_file_size_bytes(const std::string& text_name) {
  uint64_t file_size_bytes = 0;
  {
    std::ifstream is(text_name, std::ios::binary | std::ios::ate);
    assert(is.good());
    file_size_bytes = is.tellg();
    is.close();
  }
  return file_size_bytes;
}

template <typename T>
std::string entry(const T& value) {
  // Convert to string.
  std::stringstream ss;
  if (std::is_same<T, double>::value) {
    ss << std::fixed << std::setprecision(4);
  }
  ss << value;
  std::string s = ss.str();

#if defined(CSV_FORMAT)

  // Escape double quotes.
  uint64_t pos = 0;
  while ((pos = s.find('"', pos)) != std::string::npos) {
    s.replace(pos, 1, "\"\"");
    pos += 2;
  }

  // Add quotes around entry.
  return "\"" + s + "\"";

#elif defined(JSON_FORMAT)

  if (
      std::is_same<T, std::string>::value 
      || std::is_same<typename std::decay<T>::type, const char*>::value
      || std::is_same<typename std::decay<T>::type, char*>::value) {

    // Escape double quotes.
    uint64_t pos = 0;
    while ((pos = s.find('"', pos)) != std::string::npos) {
      s.replace(pos, 1, "\\\"");
      pos += 2;
    }

    // Add quotes around entry.
    return "\"" + s + "\"";

  }
  else {
    return s;
  }

#else

  uint64_t width = 25;
  if (s.size() == 0) {
    s = "NA";
  }
  std::stringstream oss;
  oss << " " << std::setw(width) << s;
  return oss.str();
  
#endif

}

void print_header() {
  static bool print_header = true;
  if (print_header) {
    print_header = false;
#if defined(CSV_FORMAT)
    uint64_t n = 0;
#define X(arg) ,#arg
    auto fields = { EXPAND_AND_REMOVE_FIRST(BENCHMARK_FIELDS) };
#undef X
    for (const auto& field : fields) {
      std::cout << (n++ == 0 ? "" : ",") << entry(field);
    }
    std::cout << std::endl;
#elif defined(JSON_FORMAT)
    std::cout << "[" << std::endl;
#else
#define X(arg) ,#arg
    auto fields = { EXPAND_AND_REMOVE_FIRST(BENCHMARK_FIELDS) };
#undef X
    for (const auto& field : fields ) {
      std::cout << entry(field);
    }
    std::cout << std::endl;
#endif
  }
}

void print_body(const std::map<std::string, std::string>& map) {
#if defined(CSV_FORMAT)
  uint64_t n = 0;
#define X(arg) ,#arg
  auto fields = { EXPAND_AND_REMOVE_FIRST(BENCHMARK_FIELDS) };
#undef X
  for (const auto& field : fields) {
    std::cout << (n++ == 0 ? "" : ",");
    auto it = map.find(entry(field));
    if (it != map.end()) {
      std::cout << it->second;
    }
    else {
      std::cout << entry("");
    }
  }
  std::cout << std::endl;
#elif defined(JSON_FORMAT)
  std::cout << "{";
  uint64_t n = 0;
#define X(arg) ,#arg
  auto fields = { EXPAND_AND_REMOVE_FIRST(BENCHMARK_FIELDS) };
#undef X
  for (const auto& field : fields) {
    auto it = map.find(entry(field));
    if (it != map.end()) {
      std::cout << (n++ == 0 ? "" : ",") << it->first << ":" << it->second;
    }
  }
  std::cout << "}," << std::endl;
#else
#define X(arg) ,#arg
  auto fields = { EXPAND_AND_REMOVE_FIRST(BENCHMARK_FIELDS) };
#undef X
  for (const auto& field : fields) {
    auto it = map.find(entry(field));
    if (it != map.end()) {
      std::cout << it->second;
    }
    else {
      std::cout << entry("");
    }
  }
  std::cout << std::endl;
#endif
}

void print_footer() {
  static bool print_footer = true;
  if (print_footer) {
    print_footer = false;
#if defined(CSV_FORMAT)
#elif defined(JSON_FORMAT)
    std::cout << "{}" << std::endl;
    std::cout << "]" << std::endl;
#else
#endif
  }
}

void print_results(const std::map<std::string, std::string>& map) {
  // Print header.
#ifndef HEADERLESS
  print_header();
#endif

  // Print body.
  print_body(map);

  // Print footer.
#ifndef FOOTERLESS
  print_footer();
#endif
}
