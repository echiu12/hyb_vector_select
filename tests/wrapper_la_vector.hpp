#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iterator>

#include <sdsl/bit_vectors.hpp>

#include <la_vector.hpp>

using namespace sdsl;

// Forward Declarations.
template<uint8_t t_b, typename t_bit_vector> class LAVectorWrapperRank;
template<uint8_t t_b, typename t_bit_vector> class LAVectorWrapperSelect;

// PopIterator
template <typename t_value = uint64_t, typename t_bv = sdsl::bit_vector>
class PopIterator {

public:

  using bit_vector_type = t_bv;
  using rank_support_type = typename bit_vector_type::rank_1_type;
  using select_support_type = typename bit_vector_type::select_1_type;
  using size_type = typename t_bv::size_type;

  using iterator_category = std::random_access_iterator_tag;
  using value_type = t_value;
  using difference_type = std::ptrdiff_t;
  using reference = value_type&;
  using pointer = value_type*;

private:

  const bit_vector_type *m_bv;
  const rank_support_type *m_rank_support;
  const select_support_type *m_select_support;
  size_type m_i;
  size_type m_pops;
  value_type m_storage;

public:

  // Constructor.
  PopIterator(
      const bit_vector_type *bv = nullptr,
      const rank_support_type *rank_support = nullptr,
      const select_support_type *select_support = nullptr,
      size_type i = 0
  ) : m_bv(bv),
      m_rank_support(rank_support),
      m_select_support(select_support),
      m_i(i),
      m_pops(rank_support->rank(bv->size())),
      m_storage(0) {
    set_i(m_i);
  }

  // TODO: copy and move constructors

  // Dereference operators.
  value_type operator*() const { return m_storage; }
  pointer operator->() const { return &m_storage; }

  // Increment operators.
  PopIterator& operator++() { set_i(m_i + 1); return *this; }
  PopIterator operator++(int) { PopIterator tmp = *this; ++(*this); return tmp; }

  // Decrement operators.
  PopIterator& operator--() { set_i(m_i - 1); return *this; }
  PopIterator operator--(int) { PopIterator tmp = *this; --(*this); return tmp; }

  // Arithmetic operators.
  PopIterator& operator+=(difference_type n) { set_i(m_i + n); return *this; }
  PopIterator& operator-=(difference_type n) { set_i(m_i - n); return *this; }
  PopIterator operator+(difference_type n) const { return PopIterator(m_bv, m_rank_support, m_select_support, m_i + n); }
  PopIterator operator-(difference_type n) const { return PopIterator(m_bv, m_rank_support, m_select_support,  m_i - n); }
  difference_type operator-(const PopIterator& other) const { return m_i - other.m_i; }

  // Comparison operators.
  bool operator==(const PopIterator& other) const { return m_i == other.m_i; }
  bool operator!=(const PopIterator& other) const { return m_i != other.m_i; }
  bool operator<(const PopIterator& other) const { return m_i < other.m_i; }
  bool operator>(const PopIterator& other) const { return m_i > other.m_i; }
  bool operator<=(const PopIterator& other) const { return m_i <= other.m_i; }
  bool operator>=(const PopIterator& other) const { return m_i >= other.m_i; }

  // Subscript operator.
  value_type operator[](difference_type i) { set_i(i); return m_storage; }

private:

  void set_i(uint64_t i) {
    m_i = i;
    if (m_bv != nullptr && m_i < m_pops) {
      m_storage = m_select_support->select(m_i + 1);
    }
  }

};

// LAVectorWrapper.
template <typename t_bv = la_vector<uint32_t>>
class LAVectorWrapper {

public:

  using size_type = uint64_t;
  using value_type = uint64_t;
  using difference_type = ptrdiff_t;
  using iterator = sdsl::random_access_const_iterator<LAVectorWrapper>;
  using rank_1_type = LAVectorWrapperRank<1, t_bv>;
  using rank_0_type = LAVectorWrapperRank<0, t_bv>;
  using select_1_type = LAVectorWrapperSelect<1, t_bv>;
  using select_0_type = LAVectorWrapperSelect<0, t_bv>;

  using bit_vector_type = t_bv;

private:

  size_type m_size;
  bit_vector_type m_bv;

public:

  // Constructors.
  LAVectorWrapper() = default;
  LAVectorWrapper(const LAVectorWrapper& wrapper) {
    m_size = wrapper.m_size;
    m_bv = bit_vector_type(wrapper.begin(), wrapper.end());
  }
  LAVectorWrapper(LAVectorWrapper&& wrapper) {
    m_size = std::move(wrapper.m_size);
    m_bv = std::move(wrapper.m_bv);
  }

  // Construct from SDSL bit_vector.
  LAVectorWrapper(const sdsl::bit_vector& sdsl_bv) {
    using rank_support_type = typename sdsl::bit_vector::rank_1_type;
    using select_support_type = typename sdsl::bit_vector::select_1_type;
    rank_support_type rank_support(&sdsl_bv);
    select_support_type select_support(&sdsl_bv);
    uint64_t popcount = rank_support.rank(sdsl_bv.size());
    PopIterator<uint64_t, sdsl::bit_vector> begin(&sdsl_bv, &rank_support, &select_support, 0);
    PopIterator<uint64_t, sdsl::bit_vector> end(&sdsl_bv, &rank_support, &select_support, popcount);
    m_size = sdsl_bv.size();
    m_bv = bit_vector_type(begin, end);
  }

  // Assignment operators.
  LAVectorWrapper& operator=(const LAVectorWrapper& wrapper) {
    m_size = wrapper.m_size;
    m_bv = bit_vector_type(wrapper.m_bv);
  }
  LAVectorWrapper& operator=(LAVectorWrapper&& wrapper) {
    if (this != &wrapper) {
      m_size = std::move(wrapper.m_size);
      m_bv = std::move(wrapper.m_bv);
    }
    return *this;
  }

  // Rank and select.
  size_type rank0(size_type i) const { return i - m_bv.rank(i); }
  size_type rank1(size_type i) const { return m_bv.rank(i); }
  size_type select0(size_type i) const {
    size_type n = m_size;
    size_type ones = m_bv.rank(n);
    size_type lb = 1;
    size_type rb = ones + 1;
    size_type r0 = 0;
    size_type pos = (size_type)-1;
    // Adapted from implementation of select on SDSL's SD vector.
    while (lb < rb) {
      size_type mid = (lb + rb) / 2;
      size_type x = select1(mid);
      size_type rank0 = x + 1 - mid;
      if (rank0 >= i) {
        rb = mid;
      }
      else {
        r0 = rank0;
        pos = x;
        lb = mid + 1;
      }
    }
    return pos + i - r0;
  }
  size_type select1(size_type i) const { return m_bv.select(i); }

  // Swap.
  void swap(LAVectorWrapper& wrapper) {
    if (this != &wrapper) {
      std::swap(m_bv, wrapper.m_bv);
    }
  }

  // Iterator.
  iterator begin() const {
    return m_bv.begin();
  }

  iterator end() const {
    return m_bv.end();
  }

  // Size (note that the original la_vector size is the popcount).
  size_type size() const {
    return m_size;
  }

  // Serialize.
  size_type serialize(std::ostream& out, structure_tree_node* v = nullptr, std::string name = "") const {
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;
    written_bytes += write_member(m_size, out, child, "size");
    written_bytes += m_bv.serialize(out, v, "bv");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
  }

  // Load.
  void load(std::istream& in) {
    read_member(m_size, in);
    m_bv.load(in);
  }

};

// LAVectorWrapperRank
template<uint8_t t_b, typename t_bit_vector>
class LAVectorWrapperRank {

public:
  using bit_vector_type = LAVectorWrapper<t_bit_vector>;
  using size_type = typename bit_vector_type::size_type;

private:
  const bit_vector_type* m_v;

public:
  explicit LAVectorWrapperRank(const bit_vector_type* v = nullptr) {
    set_vector(v);
  }

  size_type rank(size_type i) const {
    if constexpr (t_b == 0) {
      return m_v->rank0(i);
    }
    else {
      return m_v->rank1(i);
    }
  }

  size_type operator()(size_type i) const {
    return rank(i);
  }

  void set_vector(const bit_vector_type* v = nullptr) {
    m_v = v;
  }

  void swap(LAVectorWrapperRank&) {}

  size_type serialize(std::ostream& out, structure_tree_node* v = nullptr, std::string name = "") const {
    return 0;
  }

  void load(std::istream& in, const bit_vector_type* v = nullptr) {}

};

// LAVectorWrapperSelect
template<uint8_t t_b, typename t_bit_vector>
class LAVectorWrapperSelect {
public:
  using bit_vector_type = LAVectorWrapper<t_bit_vector>;
  using size_type = typename bit_vector_type::size_type;

private:
  const bit_vector_type* m_v;

public:
  explicit LAVectorWrapperSelect(const bit_vector_type* v = nullptr) {
    set_vector(v);
  }

  size_type select(size_type i) const {
    if constexpr (t_b == 0) {
      return m_v->select0(i);
    }
    else {
      return m_v->select1(i);
    }
  }

  size_type operator()(size_type i) const {
    return select(i);
  }

  void set_vector(const bit_vector_type* v = nullptr) {
    m_v = v;
  }

  void swap(LAVectorWrapperSelect&) {}

  size_type serialize(std::ostream& out, structure_tree_node* v = nullptr, std::string name = "") const {
    return 0;
  }

  void load(std::istream& in, const bit_vector_type* v = nullptr) {}

};
