#pragma once

#include <cstdint>

#include <sdsl/hyb_vector.hpp>

using namespace sdsl;

// Forward Declarations.
template <typename t_bit_vector = sdsl::hyb_vector<>> class NaiveHybVectorWrapper;
template<uint8_t t_b, typename t_bit_vector> class NaiveHybVectorWrapperRank;
template<uint8_t t_b, typename t_bit_vector> class NaiveHybVectorWrapperSelect;

// NaiveHybVectorWrapper.
template <typename t_bit_vector>
class NaiveHybVectorWrapper {

public:

  using size_type = uint64_t;
  using value_type = uint64_t;
  using difference_type = ptrdiff_t;
  using iterator = random_access_const_iterator<NaiveHybVectorWrapper>;
  using rank_1_type = NaiveHybVectorWrapperRank<1, t_bit_vector>;
  using rank_0_type = NaiveHybVectorWrapperRank<0, t_bit_vector>;
  using select_1_type = NaiveHybVectorWrapperSelect<1, t_bit_vector>;
  using select_0_type = NaiveHybVectorWrapperSelect<0, t_bit_vector>;

  using bit_vector_type = t_bit_vector;
  using rank_0_support_type = typename t_bit_vector::rank_0_type;
  using rank_1_support_type = typename t_bit_vector::rank_1_type;
  using select_0_support_type = typename t_bit_vector::select_0_type;
  using select_1_support_type = typename t_bit_vector::select_1_type;

private:

  bit_vector_type m_bv;
  rank_0_support_type m_rank_0_support;
  rank_1_support_type m_rank_1_support;
  select_0_support_type m_select_0_support;
  select_1_support_type m_select_1_support;

public:

  // Constructors.
  NaiveHybVectorWrapper() = default;
  NaiveHybVectorWrapper(const NaiveHybVectorWrapper& wrapper) {
    m_bv = bit_vector_type(wrapper.m_bv);
    m_rank_0_support = rank_0_support_type(&m_bv);
    m_rank_1_support = rank_1_support_type(&m_bv);
    m_select_0_support = select_0_support_type(&m_rank_0_support);
    m_select_1_support = select_1_support_type(&m_rank_1_support);
  }
  NaiveHybVectorWrapper(NaiveHybVectorWrapper&& wrapper) { *this = std::move(wrapper); }

  // Construct from SDSL bit_vector.
  NaiveHybVectorWrapper(const sdsl::bit_vector& sdsl_bv) {
    m_bv = bit_vector_type(sdsl_bv);
    m_rank_0_support = rank_0_support_type(&m_bv);
    m_rank_1_support = rank_1_support_type(&m_bv);
    m_select_0_support = select_0_support_type(&m_rank_0_support);
    m_select_1_support = select_1_support_type(&m_rank_1_support);
  }

  // Assignment operators.
  NaiveHybVectorWrapper& operator=(const NaiveHybVectorWrapper& wrapper) {
    m_bv = bit_vector_type(wrapper.m_bv);
    m_rank_0_support = rank_0_support_type(&m_bv);
    m_rank_1_support = rank_1_support_type(&m_bv);
    m_select_0_support = select_0_support_type(&m_rank_0_support);
    m_select_1_support = select_1_support_type(&m_rank_1_support);
  }
  NaiveHybVectorWrapper& operator=(NaiveHybVectorWrapper&& wrapper) {
    if (this != &wrapper) {
      m_bv = std::move(wrapper.m_bv);
      m_rank_0_support = rank_0_support_type(&m_bv);
      m_rank_1_support = rank_1_support_type(&m_bv);
      m_select_0_support = select_0_support_type(&m_rank_0_support);
      m_select_1_support = select_1_support_type(&m_rank_1_support);
    }
    return *this;
  }

  // Rank and select.
  size_type rank0(size_type i) const { return m_rank_0_support(i); }
  size_type rank1(size_type i) const { return m_rank_1_support(i); }
  size_type select0(size_type i) const { return m_select_0_support(i); }
  size_type select1(size_type i) const { return m_select_1_support(i); }

  // Swap.
  void swap(NaiveHybVectorWrapper& wrapper) {
    if (this != &wrapper) {
      std::swap(m_bv, wrapper.m_bv);
      sdsl::util::swap_support(m_rank_0_support, wrapper.m_rank_0_support, &m_bv, &wrapper.m_bv);
      sdsl::util::swap_support(m_rank_1_support, wrapper.m_rank_1_support, &m_bv, &wrapper.m_bv);
      m_select_0_support = select_0_support_type(&m_rank_0_support);
      m_select_1_support = select_1_support_type(&m_rank_1_support);
      wrapper.m_select_0_support = select_0_support_type(&wrapper.m_rank_0_support);
      wrapper.m_select_1_support = select_1_support_type(&wrapper.m_rank_1_support);
    }
  }

  // Iterator.
  iterator begin() const {
    return iterator(&m_bv, 0);
  }

  iterator end() const {
    return iterator(&m_bv, m_bv.size());
  }

  // Serialize.
  size_type serialize(std::ostream& out, structure_tree_node* v = nullptr, std::string name = "") const {
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;
    written_bytes += m_bv.serialize(out, v, "bv");
    written_bytes += m_rank_0_support.serialize(out, v, "rank_0_support");
    written_bytes += m_rank_1_support.serialize(out, v, "rank_1_support");
    written_bytes += m_select_0_support.serialize(out, v, "select_0_support");
    written_bytes += m_select_1_support.serialize(out, v, "select_1_support");
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
  }

  // Load.
  void load(std::istream& in) {
    m_bv.load(in);
    m_rank_0_support.load(in);
    m_rank_1_support.load(in);
    m_select_0_support.load(in);
    m_select_1_support.load(in);
  }

};

// NaiveHybVectorWrapperRank
template<uint8_t t_b, typename t_bit_vector>
class NaiveHybVectorWrapperRank {

public:
  using bit_vector_type = NaiveHybVectorWrapper<t_bit_vector>;
  using size_type = typename bit_vector_type::size_type;

private:
  const bit_vector_type* m_v;

public:
  explicit NaiveHybVectorWrapperRank(const bit_vector_type* v = nullptr) {
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

  void swap(NaiveHybVectorWrapperRank&) {}

  size_type serialize(std::ostream& out, structure_tree_node* v = nullptr, std::string name = "") const {
    return 0;
  }

  void load(std::istream& in, const bit_vector_type* v = nullptr) {}

};

// NaiveHybVectorWrapperSelect
template<uint8_t t_b, typename t_bit_vector>
class NaiveHybVectorWrapperSelect {
public:
  using bit_vector_type = NaiveHybVectorWrapper<t_bit_vector>;
  using size_type = typename bit_vector_type::size_type;

private:
  const bit_vector_type* m_v;

public:
  explicit NaiveHybVectorWrapperSelect(const bit_vector_type* v = nullptr) {
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

  void swap(NaiveHybVectorWrapperSelect&) {}

  size_type serialize(std::ostream& out, structure_tree_node* v = nullptr, std::string name = "") const {
    return 0;
  }

  void load(std::istream& in, const bit_vector_type* v = nullptr) {}

};
