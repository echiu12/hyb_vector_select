// 20260124_132500: included the wide_rank_select support.

#pragma once

#include <sdsl/bit_vectors.hpp>

#include <pasta/bit_vector/bit_vector.hpp>
#include <pasta/bit_vector/support/flat_rank_select.hpp>
#include <pasta/bit_vector/support/wide_rank_select.hpp>

using namespace sdsl;

// Forward Declarations.
template <typename t_bit_vector = pasta::BitVector, typename t_rank_select_support = pasta::FlatRankSelect<>> class PastaVectorWrapper;
template<uint8_t t_b, typename t_bit_vector, typename t_rank_select_support> class PastaVectorWrapperRank;
template<uint8_t t_b, typename t_bit_vector, typename t_rank_select_support> class PastaVectorWrapperSelect;

// PastaVectorWrapper.
template <typename t_bit_vector, typename t_rank_select_support>
class PastaVectorWrapper {

public:

  using size_type = uint64_t;
  using value_type = uint64_t;
  using difference_type = ptrdiff_t;
  using iterator = random_access_const_iterator<PastaVectorWrapper>;
  using rank_1_type = PastaVectorWrapperRank<1, t_bit_vector, t_rank_select_support>;
  using rank_0_type = PastaVectorWrapperRank<0, t_bit_vector, t_rank_select_support>;
  using select_1_type = PastaVectorWrapperSelect<1, t_bit_vector, t_rank_select_support>;
  using select_0_type = PastaVectorWrapperSelect<0, t_bit_vector, t_rank_select_support>;

  using bit_vector_type = t_bit_vector;
  using rank_select_support_type = t_rank_select_support;

private:

  bit_vector_type m_bv;
  rank_select_support_type m_support;

public:

  // Constructors.
  PastaVectorWrapper() = default;
  PastaVectorWrapper(const PastaVectorWrapper& wrapper) {
    m_bv = bit_vector_type(wrapper.m_bv);
    m_support = rank_select_support_type(m_bv);
  }
  PastaVectorWrapper(PastaVectorWrapper&& wrapper) { *this = std::move(wrapper); }

  // Construct from SDSL bit_vector.
  PastaVectorWrapper(const sdsl::bit_vector& sdsl_bv) {
    uint64_t n = sdsl_bv.size();
    bit_vector_type bv(n);
    for (uint64_t i = 0; i < n; ++i) {
      bv[i] = sdsl_bv[i];
    }
    m_bv = std::move(bv);
    m_support = rank_select_support_type(m_bv);
  }

  // Assignment operators.
  PastaVectorWrapper& operator=(const PastaVectorWrapper& wrapper) {
    m_bv = bit_vector_type(wrapper.m_bv);
    m_support = rank_select_support_type(m_bv);
  }
  PastaVectorWrapper& operator=(PastaVectorWrapper&& wrapper) {
    if (this != &wrapper) {
      m_bv = std::move(wrapper.m_bv);
      m_support = rank_select_support_type(m_bv);
    }
    return *this;
  }

  // Rank and select.
  size_type rank0(size_type i) const { return m_support.rank0(i); }
  size_type rank1(size_type i) const { return m_support.rank1(i); }
  size_type select0(size_type i) const { return m_support.select0(i); }
  size_type select1(size_type i) const { return m_support.select1(i); }

  // Swap.
  void swap(PastaVectorWrapper& wrapper) {
    if (this != &wrapper) {
      std::swap(m_bv, wrapper.m_bv);
      std::swap(m_support, wrapper.m_support);
      m_support = rank_select_support_type(m_bv);
      wrapper.m_support = rank_select_support_type(wrapper.m_bv);
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
    // Restore sdsl bitvector.
    uint64_t n = m_bv.size();
    sdsl::bit_vector sdsl_bv(n);
    for (uint64_t i = 0; i < n; ++i) {
      sdsl_bv[i] = m_bv[i];
    }

    // Write bitvector.
    structure_tree_node* child = structure_tree::add_child(v, name, util::class_name(*this));
    size_type written_bytes = 0;
    written_bytes += sdsl_bv.serialize(out, v, "bv");
    written_bytes += m_support.space_usage();
    structure_tree::add_size(child, written_bytes);
    return written_bytes;
  }

  // Load.
  void load(std::istream& in) {
    sdsl::bit_vector sdsl_bv;
    sdsl_bv.load(in);

    uint64_t n = sdsl_bv.size();
    pasta::BitVector bv(n);
    for (uint64_t i = 0; i < n; ++i) {
      bv[i] = sdsl_bv[i];
    }
    m_bv = std::move(bv);
  }

};

// PastaVectorWrapperRank
template<uint8_t t_b, typename t_bit_vector, typename t_rank_select_support>
class PastaVectorWrapperRank {

public:
  using bit_vector_type = PastaVectorWrapper<t_bit_vector, t_rank_select_support>;
  using size_type = typename bit_vector_type::size_type;

private:
  const bit_vector_type* m_v;

public:
  explicit PastaVectorWrapperRank(const bit_vector_type* v = nullptr) {
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

  void swap(PastaVectorWrapperRank&) {}

  size_type serialize(std::ostream& out, structure_tree_node* v = nullptr, std::string name = "") const {
    return 0;
  }

  void load(std::istream& in, const bit_vector_type* v = nullptr) {}

};

// PastaVectorWrapperSelect
template<uint8_t t_b, typename t_bit_vector, typename t_rank_select_support>
class PastaVectorWrapperSelect {
public:
  using bit_vector_type = PastaVectorWrapper<t_bit_vector, t_rank_select_support>;
  using size_type = typename bit_vector_type::size_type;

private:
  const bit_vector_type* m_v;

public:
  explicit PastaVectorWrapperSelect(const bit_vector_type* v = nullptr) {
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

  void swap(PastaVectorWrapperSelect&) {}

  size_type serialize(std::ostream& out, structure_tree_node* v = nullptr, std::string name = "") const {
    return 0;
  }

  void load(std::istream& in, const bit_vector_type* v = nullptr) {}

};
