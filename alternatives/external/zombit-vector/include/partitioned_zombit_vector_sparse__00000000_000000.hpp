/***
BSD 2-Clause License

Copyright (c) 2018, Adrián
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**/


//
// Created by Adrián on 17/12/2018.
//

#ifndef RUNS_VECTORS_PARTITIONED_ZOMBIT_VECTOR_SPARSE_HPP
#define RUNS_VECTORS_PARTITIONED_ZOMBIT_VECTOR_SPARSE_HPP

#include <sdsl/int_vector.hpp>
#include <sdsl/rank_support.hpp>
#include <algorithm>
#include <cstdint>
#include <sdsl/succ_support_v.hpp>
#include <sdsl/select_support.hpp>
#include <util.hpp>
#include <optimal_byte_partition.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/math.hpp>


#define ZOMBIT_SMALL_BLOCK 8

namespace runs_vectors {

    template<uint8_t t_b, class t_mixed> class rank_support_partitioned_zombit_sparse_v;
    template<uint8_t t_b, class t_mixed> class select_support_partitioned_zombit_sparse;
    class succ_support_partitioned_zombit_sparse_naive;
    template<uint8_t t_b, class t_mixed> class succ_support_partitioned_zombit_sparse;
    class partitioned_zombit_sparse_iterator;


    template <class t_mixed = sdsl::bit_vector>
    class partitioned_zombit_vector_sparse {

    public:
        typedef sdsl::bit_vector    bitmap_type;
        typedef typename sdsl::bit_vector::rank_1_type rank_bitmap_type;
        typedef sdsl::sd_vector<>   length_type;
        typedef typename sdsl::sd_vector<>::rank_1_type   rank_length_type;
        typedef typename sdsl::sd_vector<>::select_1_type select_length_type;
        typedef t_mixed  mixed_bitmap_type;
        typedef typename sdsl::bit_vector::size_type size_type;
        typedef typename sdsl::bit_vector::value_type value_type;
        typedef typename sdsl::bit_vector::difference_type difference_type;
        typedef partitioned_zombit_sparse_iterator plain_iterator_type;
        typedef rank_support_partitioned_zombit_sparse_v<1, t_mixed> rank_1_type;
        typedef rank_support_partitioned_zombit_sparse_v<0, t_mixed> rank_0_type;
        typedef select_support_partitioned_zombit_sparse<1, t_mixed> select_1_type;
        typedef select_support_partitioned_zombit_sparse<0, t_mixed> select_0_type;
        typedef succ_support_partitioned_zombit_sparse<1, t_mixed> succ_1_type;
        typedef sdsl::random_access_const_iterator<partitioned_zombit_vector_sparse> iterator;

        friend class select_support_partitioned_zombit_sparse<1, t_mixed>;
        friend class select_support_partitioned_zombit_sparse<0, t_mixed>;
        friend class rank_support_partitioned_zombit_sparse_v<1, t_mixed>;
        friend class rank_support_partitioned_zombit_sparse_v<0, t_mixed>;
        friend class succ_support_partitioned_zombit_sparse_naive;
        friend class succ_support_partitioned_zombit_sparse<1, t_mixed>;
        friend class succ_support_partitioned_zombit_sparse<0, t_mixed>;
        friend class partitioned_zombit_sparse_iterator;

    private:

        size_type m_size = 0;
        length_type m_partitions;
        rank_length_type m_rank_partitions;
        select_length_type m_select_partitions;
        bitmap_type m_full;
        rank_bitmap_type m_rank_full;
        bitmap_type m_info;
        rank_bitmap_type m_rank_info;
        length_type m_length_mixed;
        select_length_type m_select_length_mixed;
        mixed_bitmap_type m_mixed;

        void copy(const partitioned_zombit_vector_sparse &o){
            m_size = o.m_size;
            m_partitions = o.m_partitions;
            m_rank_partitions = o.m_rank_partitions;
            m_select_partitions = o.m_select_partitions;
            m_rank_partitions.set_vector(&m_partitions);
            m_select_partitions.set_vector(&m_partitions);
            m_full = o.m_full;
            m_rank_full = o.m_rank_full;
            m_rank_full.set_vector(&m_full);
            m_info = o.m_info;
            m_rank_info = o.m_rank_info;
            m_rank_info.set_vector(&m_info);
            m_length_mixed = o.m_length_mixed;
            m_select_length_mixed = o.m_select_length_mixed;
            m_select_length_mixed.set_vector(&m_length_mixed);
            m_mixed = o.m_mixed;
        }

        void get_small_blocks(const sdsl::bit_vector &c, std::vector<partition::block_t> &blocks_vec){
            const uint8_t* data = (uint8_t*) c.data();
            uint64_t blocks = sdsl::util::math::ceil_div(c.size(), ZOMBIT_SMALL_BLOCK);
            blocks_vec.resize(blocks);
            uint64_t r0 = 0, r1=0, m = 0;
            for(uint64_t i = 0; i < blocks; ++i){
                if(i < blocks-1){
                    if(data[i] == 0){
                        blocks_vec[i] = partition::block_t::run0;
                        ++r0;
                    }else if (data[i] == 255){
                        blocks_vec[i] = partition::block_t::run1;
                        ++r1;
                    }else{
                        blocks_vec[i] = partition::block_t::mixed;
                        ++m;
                    }
                }else{
                    uint8_t mask = sdsl::bits::lo_set[c.size()-i*ZOMBIT_SMALL_BLOCK];
                    uint8_t aux = data[i] & mask;
                    if(aux == 0){
                        blocks_vec[i] = partition::block_t::run0;
                        ++r0;
                    }else if (aux == mask){
                        blocks_vec[i] = partition::block_t::run1;
                        ++r1;
                    }else{
                        blocks_vec[i] = partition::block_t::mixed;
                        ++m;
                    }
                }
            }

            std::cout << std::endl << "blocks: " << blocks << " r0: " << r0 << " r1: " << r1 << " m: " << m << std::endl;
        }

        inline size_type pos_to_block(const size_type pos) const{
            return m_rank_partitions(pos/ZOMBIT_SMALL_BLOCK+1) - 1;
        }

        inline size_type begin_block(const size_type b) const{
            return m_select_partitions(b+1) * ZOMBIT_SMALL_BLOCK;
        }

        inline size_type begin_mixed_block(const size_type mb) const{
            return m_select_length_mixed(mb+1) * ZOMBIT_SMALL_BLOCK;
        }


    public:

        const bitmap_type &full = m_full;
        const bitmap_type &info = m_info;
        const mixed_bitmap_type &mixed = m_mixed;
        const select_length_type &select_partitions = m_select_partitions;
        const select_length_type &select_length_mixed = m_select_length_mixed;

        partitioned_zombit_vector_sparse() = default;



        explicit partitioned_zombit_vector_sparse(const sdsl::bit_vector &c){

            if(c.empty()) return;

            //Computing runs and sample
            m_size = c.size();
            std::vector<partition::block_t> small_blocks;
            get_small_blocks(c, small_blocks);
            partition::optimal_sparse_variable opt_var(small_blocks.begin(), small_blocks.size(), c.size());
            size_type n_blocks = opt_var.partition.size();

            //Computing full and mixed blocks
            m_full = sdsl::bit_vector(n_blocks+1, 0);
            m_info = sdsl::bit_vector(n_blocks+1, 0);

            //Computing full and mixed blocks
            size_type i_sb = 0, cnt_mixed_blocks = 0, i_f = 0;
            std::vector<std::pair<size_type, size_type>> mixed_blocks;
           // std::vector<bool> is_mixed(n_blocks, false);
            for(size_type i_p = 0; i_p < n_blocks; ++i_p) {
                size_type cnt_run0 = 0, cnt_run1 = 0, cnt_mixed = 0;
                size_type beg_sb = i_sb;
                while (i_sb < opt_var.partition[i_p]) {
                    if (small_blocks[i_sb] == partition::run0) {
                        ++cnt_run0;
                    } else if (small_blocks[i_sb] == partition::run1) {
                        ++cnt_run1;
                    } else {
                        ++cnt_mixed;
                       // is_mixed[i_sb] = true;
                    }
                    ++i_sb;
                }

                if (cnt_mixed || (cnt_run0 && cnt_run1)) { //Mixed
                    m_full[i_f++] = 0;
                    m_info[i_p] = 1;
                    cnt_mixed_blocks += opt_var.partition[i_p] - beg_sb;
                    mixed_blocks.emplace_back(beg_sb, opt_var.partition[i_p]);
                } else if (cnt_run1) { //Full Ones
                    m_full[i_f++] = 1;
                    m_info[i_p] = 1;
                } else { //Empty
                    m_info[i_p] = 0;
                }
            }
            //tricks for avoiding if statements in succ
            m_info[n_blocks] = 1;
            m_full[i_f] = 1;
            m_full.resize(i_f+1);
            //Storing info of mixed blocks
            auto aux = sdsl::bit_vector(cnt_mixed_blocks*ZOMBIT_SMALL_BLOCK);
            size_type ith_mixed = 0, start_block, end_block;
            for(const auto &mixed_block : mixed_blocks){
                start_block = mixed_block.first * ZOMBIT_SMALL_BLOCK;
                end_block = std::min(mixed_block.second * ZOMBIT_SMALL_BLOCK, c.size());
                for(size_type i = start_block; i < end_block; ++i){
                    aux[ith_mixed] = c[i];
                    ++ith_mixed;
                }
            }
            m_mixed = mixed_bitmap_type(aux);
            //Storing lengths
            {
                sdsl::bit_vector part_aux(small_blocks.size()+1, 0);
                part_aux[0]=1;
                for(size_type ith_block = 0; ith_block < n_blocks; ++ith_block){
                    part_aux[opt_var.partition[ith_block]] = 1;
                }
                //part_aux[n_blocks]=1;
                m_partitions = sdsl::sd_vector<>(part_aux);
            }


            {
                sdsl::bit_vector mix_aux(cnt_mixed_blocks+1, 0);
                mix_aux[0]=1;
                size_type last = 0;
                for(const auto &mixed_block : mixed_blocks){
                    auto len = mixed_block.second - mixed_block.first;
                    last = last + len;
                    mix_aux[last] = 1;
                }
                mix_aux[cnt_mixed_blocks]=1;
                m_length_mixed = sdsl::sd_vector<>(mix_aux);
            }

            sdsl::util::init_support(m_rank_info, &m_info);
            sdsl::util::init_support(m_rank_full, &m_full);
            sdsl::util::init_support(m_rank_partitions, &m_partitions);
            sdsl::util::init_support(m_select_partitions, &m_partitions);
            sdsl::util::init_support(m_select_length_mixed, &m_length_mixed);
        }

        //! Copy constructor
        partitioned_zombit_vector_sparse(const partitioned_zombit_vector_sparse& o)
        {
            copy(o);
        }

        //! Move constructor
        partitioned_zombit_vector_sparse(partitioned_zombit_vector_sparse&& o)
        {
            *this = std::move(o);
        }


        partitioned_zombit_vector_sparse &operator=(const partitioned_zombit_vector_sparse &o) {
            if (this != &o) {
                copy(o);
            }
            return *this;
        }
        partitioned_zombit_vector_sparse &operator=(partitioned_zombit_vector_sparse &&o) {
            if (this != &o) {
                m_size = o.m_size;
                m_partitions = std::move(o.m_partitions);
                m_rank_partitions = std::move(o.m_rank_partitions);
                m_rank_partitions.set_vector(&m_partitions);
                m_select_partitions = std::move(o.m_select_partitions);
                m_select_partitions.set_vector(&m_partitions);
                m_full = std::move(o.m_full);
                m_rank_full = std::move(o.m_rank_full);
                m_rank_full.set_vector(&m_full);
                m_rank_info = std::move(o.m_rank_info);
                m_rank_info.set_vector(&m_info);
                m_info = std::move(o.m_info);
                m_length_mixed = std::move(o.m_length_mixed);
                m_select_length_mixed = std::move(o.m_select_length_mixed);
                m_select_length_mixed.set_vector(&m_length_mixed);
                m_mixed = std::move(o.m_mixed);
            }
            return *this;
        }

        void swap(partitioned_zombit_vector_sparse &o) {
            // m_bp.swap(bp_support.m_bp); use set_vector to set the supported bit_vector
            std::swap(m_size, o.m_size);
            m_partitions.swap( o.m_partitions);
            sdsl::util::swap_support(m_rank_partitions, o.m_rank_partitions,
                                     &m_partitions, &o.m_partitions);
            sdsl::util::swap_support(m_select_partitions, o.m_select_partitions,
                                     &m_partitions, &o.m_partitions);
            m_full.swap(o.m_full);
            sdsl::util::swap_support(m_rank_full, o.m_rank_full, &m_full, &o.m_full);
            m_info.swap(o.m_info);
            sdsl::util::swap_support(m_rank_info, o.m_rank_info, &m_info, &o.m_info);
            m_length_mixed.swap(o.m_length_mixed);
            sdsl::util::swap_support(m_select_length_mixed, o.m_select_length_mixed, &m_length_mixed, &o.m_length_mixed);
            m_mixed.swap(o.m_mixed);
        }

        iterator begin() const
        {
            return iterator(this, 0);
        }

        iterator end() const
        {
            return iterator(this, size());
        }

        inline value_type access(const size_type i) const{

            auto j = pos_to_block(i);
            if(!m_info[j]) return 0;
            auto i_wo = m_rank_info(j+1) - 1;
            if(m_full[i_wo]) return 1;
            auto mix_b = i_wo - m_rank_full(i_wo + 1);
            auto delta = i - begin_block(j);
            return m_mixed[begin_mixed_block(mix_b) + delta];

        }

        inline value_type operator[](const size_type i) const{
            return access(i);
        }

        inline size_type size() const {
            return m_size;
        }

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="")const
        {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_partitions.serialize(out, child, "partitions");
            written_bytes += m_rank_partitions.serialize(out, child, "rank_partitions");
            written_bytes += m_select_partitions.serialize(out, child, "select_partitions");
            written_bytes += m_full.serialize(out, child, "full");
            written_bytes += m_rank_full.serialize(out, child, "rank_full");
            written_bytes += m_info.serialize(out, child, "info");
            written_bytes += m_rank_info.serialize(out, child, "rank_info");
            written_bytes += m_length_mixed.serialize(out, child, "length_mixed");
            written_bytes += m_select_length_mixed.serialize(out, child, "select_length_mixed");
            written_bytes += m_mixed.serialize(out, child, "mixed");
            written_bytes += sdsl::serialize(m_size, out, child, "size");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in)
        {
            m_partitions.load(in);
            m_rank_partitions.load(in, &m_partitions);
            m_select_partitions.load(in, &m_partitions);
            m_full.load(in);
            m_rank_full.load(in, &m_full);
            m_info.load(in);
            m_rank_info.load(in, &m_info);
            m_length_mixed.load(in);
            m_select_length_mixed.load(in, &m_length_mixed);
            m_mixed.load(in);
            sdsl::load(m_size, in);
        }

        /*void print(){
            std::cout << "sample=" << sample << std::endl;
            std::cout << "full={";
            for(auto i = 0; i < m_full.size(); ++i){
                std::cout << m_full[i] << ", ";
            }
            std::cout << "}" << std::endl;

            std::cout << "info={";
            for(auto i = 0; i < m_info.size(); ++i){
                std::cout << m_info[i] << ", ";
            }
            std::cout << "}" << std::endl;

            std::cout << "mixed={";
            for(auto i = 0; i < m_mixed.size(); ++i){
                std::cout << m_mixed[i] << ", ";
            }
            std::cout << "}" << std::endl;
        }*/

    };



    template<uint8_t t_b>
    struct rank_support_partitioned_zombit_sparse_trait {
        typedef partitioned_zombit_vector_sparse<>::size_type size_type;
        static size_type adjust_rank(size_type r,size_type)
        {
            return r;
        }
    };

    template<>
    struct rank_support_partitioned_zombit_sparse_trait<0> {
        typedef partitioned_zombit_vector_sparse<>::size_type size_type;
        static size_type adjust_rank(size_type r, size_type n)
        {
            return n - r;
        }
    };

    template<uint8_t t_b, class t_mixed>
    class rank_support_partitioned_zombit_sparse_v {

    public:
        typedef sdsl::bit_vector::size_type size_type;
        typedef sdsl::bit_vector::value_type value_type;
        typedef typename t_mixed::rank_1_type rank_mixed_type;
        typedef typename sdsl::bit_vector::rank_1_type rank_type;
        typedef sdsl::rank_support_trait<t_b, 1> trait_type;
    private:
        const partitioned_zombit_vector_sparse<t_mixed>* m_v = nullptr;

        rank_mixed_type m_rank_mixed;
        sdsl::sd_vector<>  m_length_o;
        sdsl::select_support_sd<1> m_select_length_o;

        void copy(const rank_support_partitioned_zombit_sparse_v& ss){
            m_v = ss.m_v;
            m_rank_mixed = ss.m_rank_mixed;
            m_length_o = ss.m_length_o;
            m_select_length_o.set_vector(&m_length_o);
            if(m_v != nullptr){
                m_rank_mixed.set_vector(&(m_v->mixed));
            }
        }
    public:


        rank_support_partitioned_zombit_sparse_v() = default;

        //! Copy constructor
        rank_support_partitioned_zombit_sparse_v(const rank_support_partitioned_zombit_sparse_v& hybrid)
        {
            copy(hybrid);
        }

        explicit rank_support_partitioned_zombit_sparse_v(const partitioned_zombit_vector_sparse<t_mixed>* v)
        {
            m_v = v;
            sdsl::util::init_support(m_rank_mixed, &(m_v->mixed));
            size_type n_partitions = m_v->info.size()-1;
            sdsl::bit_vector aux(m_v->m_partitions.size(), 0);
            aux[0]=1;
            size_type last = 0, i_f = 0;
            for(size_type i = 1; i <= n_partitions; ++i){
                if(m_v->info[i-1]){
                    if (m_v->full[i_f++]) {
                        if (i > 1) {
                            last = last + (m_v->m_select_partitions(i + 1) - m_v->m_select_partitions(i));
                        } else {
                            last = m_v->m_select_partitions(i + 1);
                        }
                        aux[last] = 1;
                    }
                }
            }
            aux.resize(last+1);
            m_length_o = sdsl::sd_vector<>(aux);
            sdsl::util::init_support(m_select_length_o, &m_length_o);
        }


        //! Returns the position of the i-th occurrence in the bit vector.
        size_type rank(size_type i)const
        {
            assert(i >= 0); assert(i <= m_v->size());
            if(i <= 0) return 0;
            auto j = m_v->pos_to_block(i-1); //block

            size_type w_o = m_v->m_rank_info(j);
            size_type o = m_v->m_rank_full(w_o); //ones blocks before j
            size_type m = w_o - o; //mixed blocks before j

            size_type r = o ? m_select_length_o(o+1)*ZOMBIT_SMALL_BLOCK : 0; //1-bits of ones blocks
            if(m_v->info[j]){
                if(m_v->full[w_o]){
                    r += (i - m_v->begin_block(j)) + m_rank_mixed(m_v->m_select_length_mixed(m+1)*ZOMBIT_SMALL_BLOCK);
                }else{
                    r += m_rank_mixed(m_v->m_select_length_mixed(m+1)*ZOMBIT_SMALL_BLOCK + (i - m_v->begin_block(j)));
                }
            }else{
                r += m_rank_mixed(m_v->m_select_length_mixed(m+1)*ZOMBIT_SMALL_BLOCK); //1-bits of mixed blocks including the current one
            }
            return rank_support_partitioned_zombit_sparse_trait<t_b>::adjust_rank(r, i);
        }


        size_type operator()(size_type i)const
        {
            return rank(i);
        }

        size_type size()const
        {
            return m_v->size();
        }

        void set_vector(const partitioned_zombit_vector_sparse<t_mixed>* v=nullptr)
        {
            m_v = v;
            if(m_v != nullptr){
                sdsl::util::init_support(m_rank_mixed, &(m_v->mixed));
            }

        }

        rank_support_partitioned_zombit_sparse_v& operator=(const rank_support_partitioned_zombit_sparse_v& ss)
        {
            if (this != &ss) {
                copy(ss);
            }
            return *this;
        }

        rank_support_partitioned_zombit_sparse_v& operator=(rank_support_partitioned_zombit_sparse_v& ss)
        {

            if (this != &ss) {
                m_v = std::move(ss.m_v);
                m_rank_mixed = std::move(ss.m_rank_mixed);
                m_length_o = std::move(ss.m_length_o);
                m_select_length_o.set_vector(&m_length_o);
                if(m_v != nullptr){
                    m_rank_mixed.set_vector(&(m_v->mixed));
                }
            }
            return *this;
        }

        void swap(rank_support_partitioned_zombit_sparse_v& o) {
            std::swap(m_v, o.m_v);
            m_rank_mixed.swap(o.m_rank_mixed);
            m_length_o.swap(o.m_length_o);
            sdsl::util::swap_support(m_select_length_o, o.m_select_length_o, &m_length_o, &o.m_length_o);
            if(m_v != nullptr){
                m_rank_mixed.set_vector(&(m_v->mixed));
            }else{
                m_rank_mixed.set_vector(nullptr);
            }
            if(o.m_v != nullptr){
                o.m_rank_mixed.set_vector(&(o.m_v->mixed));
            }else{
                o.m_rank_mixed.set_vector(nullptr);
            }


        }



        void load(std::istream& in, const partitioned_zombit_vector_sparse<t_mixed>* v=nullptr)
        {
            m_v = v;
            m_length_o.load(in);
            m_select_length_o.load(in, &m_length_o);
            if(v != nullptr){
                m_rank_mixed.load(in, &(m_v->mixed));
            }


        }

        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="")const
        {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_length_o.serialize(out, child, "length_o");
            written_bytes += m_select_length_o.serialize(out, child, "select_length_o");
            written_bytes += m_rank_mixed.serialize(out, child, "rank_mixed");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }


    };


    template<uint8_t t_b, class t_mixed = sdsl::bit_vector>
    class succ_support_partitioned_zombit_sparse {

    public:
        typedef sdsl::bit_vector::size_type size_type;
        typedef sdsl::bit_vector::value_type value_type;
        typedef typename t_mixed::succ_1_type succ_1_mixed_type;
    private:
        const partitioned_zombit_vector_sparse<t_mixed>* m_v;
        sdsl::succ_support_v<1>                   m_succ_info;
        succ_1_mixed_type                         m_succ_mixed;


        void copy(const succ_support_partitioned_zombit_sparse& ss){
            m_v = ss.m_v;
            m_succ_info = ss.m_succ_info;
            m_succ_mixed = ss.m_succ_mixed;
            if(m_v != nullptr){
                m_succ_info.set_vector(&m_v->m_info);
                m_succ_mixed.set_vector(&(m_v->m_mixed));
            }else{
                m_succ_info.set_vector(nullptr);
                m_succ_mixed.set_vector(nullptr);
            }
        }

        inline size_type next_block(size_type j, size_type w_o, size_type m) const {
            j = m_succ_info(j+1); //Next partition with 1-bits
            //if(j == m_v->m_info.size()-1) return m_v->size();
            if(m_v->m_full[w_o]) return m_v->begin_block(j);
            size_type lb = m_v->begin_mixed_block(m);
            return (m_succ_mixed(lb) - lb) + m_v->begin_block(j);
        }

        inline size_type next_block(size_type j, size_type w_o, size_type next_m, size_type lb) const {
            j = m_succ_info(j+1);
            //if(j == m_v->m_info.size()-1) return m_v->size();
            if(m_v->m_full[w_o]) return m_v->begin_block(j);
            return (next_m - lb) + m_v->begin_block(j);
        }

    public:

        //! Copy constructor
        succ_support_partitioned_zombit_sparse(const succ_support_partitioned_zombit_sparse& hybrid)
        {
            copy(hybrid);
        }

        succ_support_partitioned_zombit_sparse(const partitioned_zombit_vector_sparse<t_mixed>* v = nullptr)
        {
            m_v = v;
            if(m_v == nullptr) return;
            sdsl::util::init_support(m_succ_info, &(m_v->m_info));
            sdsl::util::init_support(m_succ_mixed, &(m_v->m_mixed));
        }



        inline size_type succ(size_type i) const{

            auto j = m_v->pos_to_block(i);
            auto w_o = m_v->m_rank_info(j+1);
            auto m = w_o - m_v->m_rank_full(w_o);

            if(m_v->m_info[j]) {
                if (m_v->m_full[w_o-1]) return i;
                size_type beg = m_v->begin_block(j);
                size_type delta = (i - beg);
                size_type lb = m_v->begin_mixed_block(m-1);
                size_type rb = m_v->begin_mixed_block(m);
                auto next_m = m_succ_mixed(lb + delta);
                if (next_m < rb) return (next_m - lb) + beg;
                return next_block(j, w_o, next_m, rb);
            }else {
                return next_block(j, w_o, m);
            }
        };

        size_type operator()(size_type i)const
        {
            return succ(i);
        }

        size_type size()const
        {
            return m_v->size();
        }

        void set_vector(const partitioned_zombit_vector_sparse<t_mixed>* v=nullptr)
        {
            m_v = v;
            if(v != nullptr){
                sdsl::util::init_support(m_succ_info, &(m_v->m_info));
                sdsl::util::init_support(m_succ_mixed, &(m_v->m_mixed));
            }
        }

        succ_support_partitioned_zombit_sparse& operator=(const succ_support_partitioned_zombit_sparse& ss)
        {
            if (this != &ss) {
                copy(ss);
            }
            return *this;
        }

        succ_support_partitioned_zombit_sparse& operator=(succ_support_partitioned_zombit_sparse&& ss)
        {
            if (this != &ss) {
                m_v = std::move(ss.m_v);
                m_succ_info = std::move(ss.m_succ_info);
                m_succ_mixed = std::move(ss.m_succ_mixed);
                if(m_v != nullptr){
                    m_succ_info.set_vector(&(m_v->m_info));
                    m_succ_mixed.set_vector(&(m_v->m_mixed));
                }else{
                    m_succ_info.set_vector(nullptr);
                    m_succ_mixed.set_vector(nullptr);
                }
            }
            return *this;
        }

        void swap(succ_support_partitioned_zombit_sparse& ss) {
            if (this != &ss) {
                std::swap(m_v, ss.m_v);
                m_succ_info.swap(ss.m_succ_info);
                m_succ_mixed.swap(ss.m_succ_mixed);
                if(m_v != nullptr){
                    m_succ_info.set_vector(&(m_v->m_info));
                    m_succ_mixed.set_vector(&(m_v->m_mixed));
                }else{
                    m_succ_info.set_vector(nullptr);
                    m_succ_mixed.set_vector(nullptr);
                }
                if(ss.m_v != nullptr){
                    ss.m_succ_info.set_vector(&(ss.m_v->m_info));
                    ss.m_succ_mixed.set_vector(&(ss.m_v->m_mixed));
                }else{
                    ss.m_succ_info.set_vector(nullptr);
                    ss.m_succ_mixed.set_vector(nullptr);
                }
            }
        }

        void load(std::istream& in, const partitioned_zombit_vector_sparse<t_mixed>* v=nullptr)
        {
            m_v = v;
            if(m_v != nullptr){
                m_succ_info.load(in, &(m_v->m_info));
                m_succ_mixed.load(in, &(m_v->m_mixed));
            }
        }

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="")const
        {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += m_succ_info.serialize(out, child, "succ_info");
            written_bytes += m_succ_mixed.serialize(out, child, "succ_mixed");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }


    };


    class partitioned_zombit_sparse_iterator {

    public:
        typedef sdsl::bit_vector::size_type size_type;
        typedef sdsl::bit_vector::value_type value_type;
    private:
        const partitioned_zombit_vector_sparse<sdsl::bit_vector>* m_v;
        size_type i_info;
        size_type i_full;
        size_type i_mixed;

        size_type b_mixed, b_block;
        size_type e_mixed, e_block;
        bool is_full;
        value_type answer;

        void copy(const partitioned_zombit_sparse_iterator& other){
            i_info = other.i_info;
            i_full = other.i_full;
            b_mixed = other.b_mixed;
            e_mixed = other.e_mixed;
            b_block = other.b_block;
            e_block = other.e_block;
            is_full = other.is_full;
            answer = other.answer;
            m_v = other.m_v;
        }

    public:
        explicit partitioned_zombit_sparse_iterator(partitioned_zombit_vector_sparse<sdsl::bit_vector>* v){
            m_v = v;
            i_info = i_full = i_mixed = answer = -1ULL;
            e_mixed = e_block = 0;

        }

        partitioned_zombit_sparse_iterator() = default;

        //! Copy constructor
        partitioned_zombit_sparse_iterator(const partitioned_zombit_sparse_iterator& other)
        {
            copy(other);
        }

        value_type operator*() const{
            return answer;
        }

        bool next(){
            if(is_full and answer < e_block-1){
                ++answer;
                return answer < m_v->size();
            }
            size_type n_m;
            if(!is_full and answer < e_block-1){
                n_m = sdsl::bits_more::next_limit(m_v->mixed.data(),
                                                  b_mixed + (answer - b_block) + 1,
                                                  e_mixed);
                if(n_m < e_mixed){
                    answer = b_block + (n_m - b_mixed);
                    return answer < m_v->size();
                }
            }
            i_info = sdsl::bits_more::next_limit(m_v->m_info.data(), i_info+1, m_v->m_info.size());
            if(i_info >= m_v->m_info.size()-1) return false;
            is_full = m_v->m_full[++i_full];
            b_block = m_v->begin_block(i_info);
            e_block = m_v->begin_block(i_info+1);
            if(is_full){
                answer = b_block;
            }else{
                ++i_mixed;
                b_mixed = e_mixed;//m_v->begin_mixed_block(i_full);
                e_mixed = m_v->begin_mixed_block(i_mixed+1);
                n_m = sdsl::bits_more::next_limit(m_v->mixed.data(),
                                                  b_mixed,
                                                  e_mixed);
                answer = b_block + (n_m - b_mixed);
            }
            return answer < m_v->size();
        }

        bool next(size_type lb){ //with lower bound
            if(lb >= m_v->size()) return false;
            if(b_block <= lb && lb < e_block){
                size_type n_m;
                if(is_full){
                    answer = lb;
                    return true;
                }else{
                    n_m = sdsl::bits_more::next_limit(m_v->mixed.data(),
                                                      b_mixed + (lb - b_block),
                                                      e_mixed);
                    if(n_m < e_mixed){
                        answer = b_block + (n_m - b_mixed);
                        //std::cout << 2 << std::endl;
                        return true;
                    }
                }
                i_info = sdsl::bits_more::next_limit(m_v->m_info.data(), i_info+1, m_v->m_info.size());
                if(i_info >= m_v->m_info.size()-1) return false;
                is_full = m_v->m_full[++i_full];
                b_block = m_v->begin_block(i_info);
                e_block = m_v->begin_block(i_info+1);
                if(is_full){
                    answer = b_block;
                }else{
                    ++i_mixed;
                    b_mixed = e_mixed;//m_v->begin_mixed_block(i_full);
                    e_mixed = m_v->begin_mixed_block(i_mixed+1);
                    n_m = sdsl::bits_more::next_limit(m_v->mixed.data(),
                                                      b_mixed,
                                                      e_mixed);
                    answer = b_block + (n_m - b_mixed);
                }
            }else{
                //Compute i_info of lb
                size_type n_m;

                auto i = m_v->pos_to_block(lb);
                i_info = sdsl::bits_more::next_limit(m_v->m_info.data(), i, m_v->m_info.size());
                if(i_info >= m_v->m_info.size()) return false;
                i_full = m_v->m_rank_info(i_info);
                is_full = m_v->m_full[i_full];
                b_block = m_v->begin_block(i_info);
                e_block = m_v->begin_block(i_info+1);
                if(is_full){
                    //std::cout << 3 << std::endl;
                    answer = (i == i_info) ? lb : b_block;
                }else{
                    i_mixed = i_full - m_v->m_rank_full(i_full+1);
                    b_mixed = m_v->begin_mixed_block(i_mixed);
                    e_mixed = m_v->begin_mixed_block(i_mixed+1);
                    n_m = sdsl::bits_more::next_limit(m_v->mixed.data(),
                                                      (i == i_info) ?  b_mixed + (lb - b_block) : b_mixed,
                                                      e_mixed);
                    answer = b_block + (n_m - b_mixed);
                    //  std::cout << 4 << std::endl;
                }
            }
            return answer < m_v->size();
        }

        partitioned_zombit_sparse_iterator& operator=(const partitioned_zombit_sparse_iterator& ss)
        {
            if (this != &ss) {
                copy(ss);
            }
            return *this;
        }

        partitioned_zombit_sparse_iterator& operator=(partitioned_zombit_sparse_iterator& other)
        {

            if (this != &other) {
                i_info = other.i_info;
                b_mixed = other.b_mixed;
                e_mixed = other.e_mixed;
                b_block = other.b_block;
                e_block = other.e_block;
                is_full = other.is_full;
                answer = other.answer;
                m_v = other.m_v;
            }
            return *this;
        }
    };

    //Only with t_b=1
    //template<uint8_t t_b, class t_mixed = sdsl::bit_vector>
    class succ_support_partitioned_zombit_sparse_naive {

    public:
        typedef sdsl::bit_vector::size_type size_type;
        typedef sdsl::bit_vector::value_type value_type;
    private:
        const partitioned_zombit_vector_sparse<sdsl::bit_vector>* m_v;


        void copy(const succ_support_partitioned_zombit_sparse_naive& ss){
            m_v = ss.m_v;
        }

        inline size_type next_block(size_type j, size_type w_o, size_type m) const {
            j = sdsl::bits_more::next_limit(m_v->m_info.data(), j+1, m_v->m_info.size()); //Next partition with 1-bits
            //if(j == m_v->m_info.size()-1) return m_v->size();
            if(m_v->m_full[w_o]) return m_v->begin_block(j);
            size_type lb = m_v->begin_mixed_block(m);
            return (sdsl::bits_more::next_limit(m_v->mixed.data(), lb, m_v->mixed.size()) - lb) + m_v->begin_block(j);
        }

        inline size_type next_block(size_type j, size_type w_o, size_type next_m, size_type lb) const {
            j = sdsl::bits_more::next_limit(m_v->m_info.data(), j+1, m_v->m_info.size());
            //if(j == m_v->m_info.size()-1) return m_v->size();
            if(m_v->m_full[w_o]) return m_v->begin_block(j);
            return (next_m - lb) + m_v->begin_block(j);
        }

   public:

        //! Copy constructor
        succ_support_partitioned_zombit_sparse_naive(const succ_support_partitioned_zombit_sparse_naive& hybrid)
        {
            copy(hybrid);
        }

        succ_support_partitioned_zombit_sparse_naive(const partitioned_zombit_vector_sparse<sdsl::bit_vector>* v = nullptr)
        {
            m_v = v;
            if(m_v == nullptr) return;
        }


        inline size_type succ(size_type i) const{

            auto j = m_v->pos_to_block(i);
            auto w_o = m_v->m_rank_info(j+1);
            auto m = w_o - m_v->m_rank_full(w_o);

            if(m_v->m_info[j]) {
                if (m_v->m_full[w_o-1]) return i;
                size_type beg = m_v->begin_block(j);
                size_type delta = (i - beg);
                size_type lb = m_v->begin_mixed_block(m-1);
                size_type rb = m_v->begin_mixed_block(m);
                auto next_m = sdsl::bits_more::next_limit(m_v->mixed.data(), lb + delta, m_v->mixed.size());
                if (next_m < rb) return (next_m - lb) + beg;
                return next_block(j, w_o, next_m, rb);
            }else {
                return next_block(j, w_o, m);
            }
        };

        size_type operator()(size_type i)const
        {
            return succ(i);
        }

        size_type size()const
        {
            return m_v->size();
        }

        void set_vector(const partitioned_zombit_vector_sparse<sdsl::bit_vector>* v=nullptr)
        {
            m_v = v;
        }

        succ_support_partitioned_zombit_sparse_naive& operator=(const succ_support_partitioned_zombit_sparse_naive& ss)
        {
            if (this != &ss) {
                copy(ss);
            }
            return *this;
        }

        succ_support_partitioned_zombit_sparse_naive& operator=(succ_support_partitioned_zombit_sparse_naive&& ss)
        {
            if (this != &ss) {
                m_v = std::move(ss.m_v);
            }
            return *this;
        }

        void swap(succ_support_partitioned_zombit_sparse_naive& ss) {
            if (this != &ss) {
                std::swap(m_v, ss.m_v);
            }
        }

        void load(std::istream& in, const partitioned_zombit_vector_sparse<sdsl::bit_vector>* v=nullptr)
        {
            m_v = v;
        }

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="")const
        {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            sdsl::structure_tree::add_size(child, 0);
            return 0;
        }


    };

    template<uint8_t t_b, class t_mixed>
    class select_support_partitioned_zombit_sparse
    {
    public:
        typedef partitioned_zombit_vector_sparse<t_mixed> bit_vector_type;
        typedef rank_support_partitioned_zombit_sparse_v<t_b, t_mixed> rank_zombit_type;
        typedef typename bit_vector_type::size_type size_type;
        enum { bit_pat = t_b };
        enum { bit_pat_len = (uint8_t)1 };
    private:
        //const bit_vector_type* m_v;
        const rank_zombit_type* m_rank;

    public:
        //! Standard constructor

        select_support_partitioned_zombit_sparse() = default;


        explicit select_support_partitioned_zombit_sparse(const partitioned_zombit_vector_sparse<t_mixed>* bv)
        {
        }

        explicit select_support_partitioned_zombit_sparse(const rank_zombit_type* rank)
        {
            m_rank = rank;
        }

        //! Answers select queries
        size_type select(size_type i) const
        {
            //fprintf(stderr, "\nzombit_vector: select queries are not currently supported\n");
            //std::exit(EXIT_FAILURE);

            uint64_t l = 0, r = m_rank->size()-1;
            uint64_t mid, cnt;
            while(l < r){
                mid = (l + r) >> 1;
                cnt = m_rank->rank(mid+1);
                if(cnt < i){
                    l = mid + 1;
                }else{
                    r = mid;
                }
            }
            return l;
        }

        //! Shorthand for select(i)
        size_type operator()(size_type i) const
        {
            return select(i);
        }

        //! Return the size of the original vector
        size_type size() const
        {
            return m_rank->size();
        }

        void set_rank(const rank_zombit_type* rank){
            m_rank = rank;
        }

        //! Assignment operator
        select_support_partitioned_zombit_sparse& operator=(const select_support_partitioned_zombit_sparse& rs)
        {
            if (this != &rs) {
                m_rank = rs.m_rank;
            }
            return *this;
        }

        //! Swap method
        void swap(select_support_partitioned_zombit_sparse&) {}

        //! Load the data structure from a stream and set the supported vector
        void load(std::istream&, const rank_zombit_type* rank)
        {
            m_rank = rank;
        }

        //! Load the data structure from a stream and set the supported vector
        void load(std::istream&, const partitioned_zombit_vector_sparse<t_mixed>* bv)
        {
        }

        void set_vector(const partitioned_zombit_vector_sparse<t_mixed>* bv){

        }

        //! Serializes the data structure into a stream
        size_type serialize(std::ostream&, sdsl::structure_tree_node* v = nullptr, std::string name = "") const
        {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            sdsl::structure_tree::add_size(child, 0);
            return 0;
        }
    };

    //! Select support for the hyb_vector class
/*!
 * \tparam t_b            The bit pattern of size one. (so `0` or `1`)
 * \tparam k_sblock_rate  Superblock rate (number of blocks inside superblock)
 * TODO: implement select queries, currently this is dummy class.
 */
  /*  template<uint8_t t_b, class t_mixed>
    class select_support_zombit_v4
    {
    public:
        typedef partitioned_zombit_vector<t_mixed> bit_vector_type;
        typedef typename bit_vector_type::size_type size_type;
        enum { bit_pat = t_b };
        enum { bit_pat_len = (uint8_t)1 };
    private:
        const bit_vector_type* m_v;

    public:
        //! Standard constructor
        explicit select_support_zombit_v4(const bit_vector_type* v = nullptr)
        {
            set_vector(v);
        }

        //! Answers select queries
        size_type select(size_type) const
        {
            fprintf(stderr, "\nzombit_vector: select queries are not currently supported\n");
            std::exit(EXIT_FAILURE);
        }

        //! Shorthand for select(i)
        const size_type operator()(size_type i) const
        {
            return select(i);
        }

        //! Return the size of the original vector
        const size_type size() const
        {
            return m_v->size();
        }

        //! Set the supported vector
        void set_vector(const bit_vector_type* v = nullptr)
        {
            m_v = v;
        }

        //! Assignment operator
        select_support_zombit_v4& operator=(const select_support_zombit_v4& rs)
        {
            if (this != &rs) {
                set_vector(rs.m_v);
            }
            return *this;
        }

        //! Swap method
        void swap(select_support_zombit_v4&) {}

        //! Load the data structure from a stream and set the supported vector
        void load(std::istream&, const bit_vector_type* v = nullptr)
        {
            set_vector(v);
        }

        //! Serializes the data structure into a stream
        size_type serialize(std::ostream&, sdsl::structure_tree_node* v = nullptr, std::string name = "") const
        {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            sdsl::structure_tree::add_size(child, 0);
            return 0;
        }
    };*/

}

#endif //RCT_ZOMBIT_VECTOR_HPP
