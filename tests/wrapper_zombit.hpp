#pragma once

#include <cstdint>

#include <sdsl/bit_vectors.hpp>

#include <zombit_vector_v3.hpp>
#include <zombit_vector_v4.hpp>
#include <partitioned_zombit_vector.hpp>
#include <partitioned_zombit_vector_sparse.hpp>
#include <rec_partitioned_zombit_vector.hpp>

using namespace sdsl;

// Forward Declarations.
template <typename t_bit_vector = runs_vectors::zombit_vector_v4<>> class ZombitVectorWrapper;
template<uint8_t t_b, typename t_bit_vector> class ZombitVectorWrapperRank;
template<uint8_t t_b, typename t_bit_vector> class ZombitVectorWrapperSelect;

// ZombitVectorWrapper.
template <typename t_bit_vector>
class ZombitVectorWrapper {

public:

  using size_type = uint64_t;
  using value_type = uint64_t;
  using difference_type = ptrdiff_t;
  using iterator = random_access_const_iterator<ZombitVectorWrapper>;
  using rank_1_type = ZombitVectorWrapperRank<1, t_bit_vector>;
  using rank_0_type = ZombitVectorWrapperRank<0, t_bit_vector>;
  using select_1_type = ZombitVectorWrapperSelect<1, t_bit_vector>;
  using select_0_type = ZombitVectorWrapperSelect<0, t_bit_vector>;

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
  ZombitVectorWrapper() = default;
  ZombitVectorWrapper(const ZombitVectorWrapper& wrapper) {
    m_bv = bit_vector_type(wrapper.m_bv);
    m_rank_0_support = rank_0_support_type(&m_bv);
    m_rank_1_support = rank_1_support_type(&m_bv);
    m_select_0_support = select_0_support_type(&m_rank_0_support);
    m_select_1_support = select_1_support_type(&m_rank_1_support);
  }
  ZombitVectorWrapper(ZombitVectorWrapper&& wrapper) { *this = std::move(wrapper); }

  // Construct from SDSL bit_vector.
  ZombitVectorWrapper(const sdsl::bit_vector& sdsl_bv) {
    m_bv = bit_vector_type(sdsl_bv);
    m_rank_0_support = rank_0_support_type(&m_bv);
    m_rank_1_support = rank_1_support_type(&m_bv);
    m_select_0_support = select_0_support_type(&m_rank_0_support);
    m_select_1_support = select_1_support_type(&m_rank_1_support);
  }

  // Assignment operators.
  ZombitVectorWrapper& operator=(const ZombitVectorWrapper& wrapper) {
    m_bv = bit_vector_type(wrapper.m_bv);
    m_rank_0_support = rank_0_support_type(&m_bv);
    m_rank_1_support = rank_1_support_type(&m_bv);
    m_select_0_support = select_0_support_type(&m_rank_0_support);
    m_select_1_support = select_1_support_type(&m_rank_1_support);
  }
  ZombitVectorWrapper& operator=(ZombitVectorWrapper&& wrapper) {
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
  void swap(ZombitVectorWrapper& wrapper) {
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

// ZombitVectorWrapperRank
template<uint8_t t_b, typename t_bit_vector>
class ZombitVectorWrapperRank {

public:
  using bit_vector_type = ZombitVectorWrapper<t_bit_vector>;
  using size_type = typename bit_vector_type::size_type;

private:
  const bit_vector_type* m_v;

public:
  explicit ZombitVectorWrapperRank(const bit_vector_type* v = nullptr) {
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

  void swap(ZombitVectorWrapperRank&) {}

  size_type serialize(std::ostream& out, structure_tree_node* v = nullptr, std::string name = "") const {
    return 0;
  }

  void load(std::istream& in, const bit_vector_type* v = nullptr) {}

};

// ZombitVectorWrapperSelect
template<uint8_t t_b, typename t_bit_vector>
class ZombitVectorWrapperSelect {
public:
  using bit_vector_type = ZombitVectorWrapper<t_bit_vector>;
  using size_type = typename bit_vector_type::size_type;

private:
  const bit_vector_type* m_v;

public:
  explicit ZombitVectorWrapperSelect(const bit_vector_type* v = nullptr) {
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

  void swap(ZombitVectorWrapperSelect&) {}

  size_type serialize(std::ostream& out, structure_tree_node* v = nullptr, std::string name = "") const {
    return 0;
  }

  void load(std::istream& in, const bit_vector_type* v = nullptr) {}

};
