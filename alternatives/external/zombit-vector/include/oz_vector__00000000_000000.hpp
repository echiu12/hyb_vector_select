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
// Created by Adrián on 14/12/2018.
//

#ifndef RUNS_VECTORS_OZ_VECTOR_HPP
#define RUNS_VECTORS_OZ_VECTOR_HPP


#include "sdsl/sd_vector.hpp"
#include "sdsl/int_vector.hpp"

namespace runs_vectors {

    template<uint8_t t_b> class rank_support_oz;
    template<uint8_t t_b> class select_support_oz;
    template<uint8_t t_b> class succ_support_oz;

    template <class t_bit_vector_o = sdsl::sd_vector<>,
            class t_bit_vector_z = sdsl::sd_vector<>,
            class t_support_rank_o = sdsl::rank_support_sd<1>,
            class t_support_rank_z = sdsl::rank_support_sd<1>,
            class t_support_select_o = sdsl::select_support_sd<1>,
            class t_support_select_z = sdsl::select_support_sd<1>>

    class oz_vector {

    public:
        typedef sdsl::bit_vector::size_type             size_type;
        typedef size_type                               value_type;
        typedef sdsl::bit_vector::difference_type       difference_type;
        typedef rank_support_oz<1> rank_1_type;
        typedef rank_support_oz<0> rank_0_type;
        typedef select_support_oz<1> select_1_type;
        typedef select_support_oz<0> select_0_type;
        typedef succ_support_oz<1> succ_1_type;
        typedef sdsl::random_access_const_iterator<oz_vector> iterator;


    private:
        t_bit_vector_o m_bit_vector_o;
        t_bit_vector_z m_bit_vector_z;
        t_support_select_o m_support_select_o;
        t_support_select_z m_support_select_z;
        t_support_rank_o m_support_rank_o;
        t_support_rank_z m_support_rank_z;
        value_type m_z_1;
        size_type m_size;
        size_type m_ones;
        size_type m_zeroes;
        size_type m_right_bsearch;

        void copy(const oz_vector<>& v)
        {
            m_bit_vector_o   = v.m_bit_vector_o;
            m_bit_vector_z  = v.m_bit_vector_z;
            m_support_select_o = v.m_support_select_o;
            m_support_select_o.set_vector(&m_bit_vector_o);
            m_support_select_z = v.m_support_select_z;
            m_support_select_z.set_vector(&m_bit_vector_z);
            m_support_rank_o = v.m_support_rank_o;
            m_support_rank_o.set_vector(&m_bit_vector_o);
            m_support_rank_z = v.m_support_rank_z;
            m_support_rank_z.set_vector(&m_bit_vector_z);
            m_z_1 = v.m_z_1;
            m_size = v.m_size;
            m_ones = v.m_ones;
            m_zeroes = v.m_zeroes;
            m_right_bsearch = v.m_right_bsearch;
        }

    public:
        const t_support_rank_o& support_rank_o = m_support_rank_o;
        const t_support_rank_z& support_rank_z = m_support_rank_z;
        const t_support_select_o& support_select_o = m_support_select_o;
        const t_support_select_z& support_select_z = m_support_select_z;
        const value_type& z_1 = m_z_1;
        const size_type& ones = m_ones;
        const size_type& zeroes = m_zeroes;
        const size_type& right_bsearch = m_right_bsearch;

        oz_vector() { };
        oz_vector(const  oz_vector& oz){
            copy(oz);
        }
        oz_vector(oz_vector&& oz){
            *this = std::move(oz);
        }

        oz_vector& operator=(const oz_vector& v)
        {
            if (this != &v) {
                copy(v);
            }
            return *this;
        }

        oz_vector& operator=(oz_vector&& v)
        {
            if (this != &v) {
                m_bit_vector_o   = std::move(v.m_bit_vector_o);
                m_bit_vector_z  = std::move(v.m_bit_vector_z);
                m_support_select_o= std::move(v.m_support_select_o);
                m_support_select_o.set_vector(&m_bit_vector_o);
                m_support_select_z = std::move(v.m_support_select_z);
                m_support_select_z.set_vector(&m_bit_vector_z);
                m_support_rank_o = std::move(v.m_support_rank_o);
                m_support_rank_o.set_vector(&m_bit_vector_o);
                m_support_rank_z = std::move(v.m_support_rank_z);
                m_support_rank_z.set_vector(&m_bit_vector_z);
                m_z_1 = std::move(v.m_z_1);
                m_size = std::move(v.m_size);
                m_ones = std::move(v.m_ones);
                m_zeroes =std::move( v.m_zeroes);
                m_right_bsearch = std::move(v.m_right_bsearch);
            }
            return *this;
        }



        oz_vector(const sdsl::bit_vector &bit_vector){

            m_ones = sdsl::util::cnt_one_bits(bit_vector);
            m_zeroes = bit_vector.size() - m_ones;
            m_ones++;
            m_zeroes++;
            sdsl::bit_vector bits_z(m_zeroes, 0);
            sdsl::bit_vector bits_o(m_ones, 0);
            m_size = bit_vector.size();
            m_z_1 = bit_vector[0] ? 0 : 1;
            size_type old_value = m_z_1;

            size_type p_o = 0;
            size_type p_z = 0;
            for(size_type i = 0; i < bit_vector.size(); i++){

                value_type value = bit_vector[i];
                if(value){
                    if(old_value != value){
                        bits_o[p_o]=1;
                    }
                    p_o++;
                }else{
                    if(old_value != value) {
                        bits_z[p_z] = 1;
                    }
                    p_z++;
                }
                old_value = value;

            }
            bits_o[p_o]=1;
            bits_z[p_z]=1;


            t_bit_vector_z bit_vector_z(bits_z);
            t_bit_vector_o bit_vector_o(bits_o);

            sdsl::util::assign(m_bit_vector_o, bit_vector_o);
            sdsl::util::assign(m_bit_vector_z, bit_vector_z);
            sdsl::util::init_support(m_support_select_o, &m_bit_vector_o);
            sdsl::util::init_support(m_support_select_z, &m_bit_vector_z);
            sdsl::util::init_support(m_support_rank_o, &m_bit_vector_o);
            sdsl::util::init_support(m_support_rank_z, &m_bit_vector_z);

            m_right_bsearch = std::max(m_support_rank_o(m_ones-1), m_support_rank_z(m_zeroes-1));// -1 bit extra

        }

        //pre: Last one must be added to itr
        template<class t_itr>
        oz_vector(const t_itr begin_0, const t_itr end_0, const t_itr begin_1, const t_itr end_1, const size_type z_1){

            if (begin_0 == end_0 || begin_1 == end_1) {
                return;
            }

            if (!is_sorted(begin_0,end_0) || !is_sorted(begin_1, end_1)) {
                throw std::runtime_error("oz_vector: source lists are not sorted.");
            }
            m_z_1 = z_1;

            t_bit_vector_z bit_vector_z(begin_0, end_0);
            t_bit_vector_o bit_vector_o(begin_1, end_1);
            m_ones = bit_vector_o.size();
            m_zeroes = bit_vector_z.size();
            m_size = bit_vector_o.size() + bit_vector_z.size()-2;

            sdsl::util::assign(m_bit_vector_o, bit_vector_o);
            sdsl::util::assign(m_bit_vector_z, bit_vector_z);
            sdsl::util::init_support(m_support_select_o, &m_bit_vector_o);
            sdsl::util::init_support(m_support_select_z, &m_bit_vector_z);
            sdsl::util::init_support(m_support_rank_o, &m_bit_vector_o);
            sdsl::util::init_support(m_support_rank_z, &m_bit_vector_z);

            m_right_bsearch = std::max(m_support_rank_o(m_ones-1), m_support_rank_z(m_zeroes-1));// -1 bit extra
        }

        size_type binary_search(size_type value)const{
            size_type l , r, m;
            l = 1; r = this->right_bsearch+1;
            size_type sum = 0;
            while (l<r) {
                m=(l+r)/2;
                sum = m_support_select_z(m) + m_support_select_o(m);
                if (sum < value)
                    l=m+1;
                else if (sum == value)
                    r=l=m+1;
                else
                    r=m;
            }
            return l-1;
        }


        value_type operator[](size_type i)const{
            size_type j = binary_search(i);
            if(m_z_1){
                if(i < m_support_select_z(j+1) + m_support_select_o(j)){
                    return 0;
                }else{
                    return 1;
                }
            }else{
                if(i < m_support_select_o(j+1) + m_support_select_z(j)){
                    return 1;
                }else{
                    return 0;
                }
            }

        }

        sdsl::bit_vector get_bit_vector() const{
            if(m_z_1){
                //starts with 0s
                sdsl::bit_vector b(m_size, 0);
                size_type i = 2;
                size_type p = support_select_z(i);
                size_type old_o = 0;
                while(p < m_size){
                    size_type select_o = support_select_o(i);
                    for(size_type j = p; j < p + select_o - old_o; j++){
                        b[j] = 1;
                    }
                    p += select_o - old_o;
                    if(p < m_size){
                        i++;
                        p = support_select_z(i) + select_o;
                        old_o = select_o;
                    }
                }
                return b;
            }else{
                //starts with 1s
                sdsl::bit_vector b(m_size, 1);
                size_type i = 2;
                size_type p = support_select_o(i);
                size_type old_z = 0;
                while(p < m_size){
                    size_type select_z = support_select_z(i);
                    for(size_type j = p; j < p + select_z - old_z; j++){
                        b[j] = 0;
                    }
                    p += select_z - old_z;
                    if(p < m_size){
                        i++;
                        p = support_select_o(i) + select_z;
                        old_z = select_z;
                    }
                }
                return b;

            }
        }

        sdsl::bit_vector get_bit_vector(const size_type l, const size_type r) const{
            if(l == r){
                return sdsl::bit_vector(8, 0);
            }
            size_type j = binary_search(l);
            size_type size = r-l;
            if(m_z_1){
                size_type select_z = m_support_select_z(j+1);
                size_type select_o = m_support_select_o(j);

                if(l < select_z + select_o){ //zeroes zone
                    sdsl::bit_vector bit_vector(size, 0); //init bit_vector to 0s
                    size_type p_bit_vector =  select_z + select_o - l; //position bit_vector
                    while(p_bit_vector < size){
                        j++;
                        select_o = m_support_select_o(j); //next 1s-run
                        while(p_bit_vector < size && l + p_bit_vector < select_z+select_o){
                            //set positions to 0
                            bit_vector[p_bit_vector] = 1;
                            p_bit_vector++;
                        }
                        if(p_bit_vector< size){
                            size_type  next_select_z = m_support_select_z(j+1); //next 0s-run
                            p_bit_vector += next_select_z - select_z; //skip |0s-run| positions
                            select_z = next_select_z;
                        }
                    }
                    return bit_vector;
                }else{
                    //l is into a run of ones
                    sdsl::bit_vector bit_vector(size, 1);
                    select_o = m_support_select_o(j+1);
                    size_type  p_bit_vector = select_o + select_z - l;
                    while(p_bit_vector < size){
                        j++;
                        select_z = m_support_select_z(j+1); //next 0s-run
                        while(p_bit_vector < size && l + p_bit_vector < select_o + select_z){
                            bit_vector[p_bit_vector] = 0;
                            p_bit_vector++;
                        }
                        if(p_bit_vector < size){
                            size_type  next_select_o = m_support_select_o(j+1); //next 1s-run
                            p_bit_vector += next_select_o - select_o; //skip |1s-run| positions
                            select_o = next_select_o;
                        }
                    }
                    return bit_vector;
                }
            }else{
                size_type select_o = m_support_select_o(j+1);
                size_type select_z = m_support_select_z(j);
                if(l < select_o + select_z){ //l is into a run of ones
                    sdsl::bit_vector bit_vector(size, 1); //init bit_vector to 1s
                    size_type p_bit_vector =  select_z + select_o - l; //position bit_vector
                    while(p_bit_vector < size){
                        j++;
                        select_z = m_support_select_z(j); //next 0s-run
                        while(p_bit_vector < size && l + p_bit_vector < select_z+select_o){
                            //set positions to 0
                            bit_vector[p_bit_vector] = 0;
                            p_bit_vector++;
                        }
                        if(p_bit_vector< size){
                            size_type  next_select_o = m_support_select_o(j+1); //next 1s-run
                            p_bit_vector += next_select_o - select_o; //skip |1s-run| positions
                            select_o = next_select_o;
                        }
                    }
                    return bit_vector;
                }else{
                    //l is into a run of zeroes
                    sdsl::bit_vector bit_vector(size, 0);
                    select_z = m_support_select_z(j+1);
                    size_type  p_bit_vector = select_o + select_z - l;
                    while(p_bit_vector < size){
                        j++;
                        select_o = m_support_select_o(j+1); //next 1s-run
                        while(p_bit_vector < size && l + p_bit_vector < select_z + select_o){
                            bit_vector[p_bit_vector] = 1;
                            p_bit_vector++;
                        }
                        if(p_bit_vector < size){
                            size_type  next_select_z = m_support_select_z(j+1); //next 0s-run
                            p_bit_vector += next_select_z - select_z; //skip |0s-run| positions
                            select_z = next_select_z;
                        }
                    }
                    return bit_vector;
                }
            }
        }

        void print() const{
            for(uint i = 0; i < this->m_bit_vector_o.size(); i++){
                std::cout << this->m_bit_vector_o[i] << ", ";
            }
            std::cout << std::endl;

            for(uint i = 0; i < this->m_bit_vector_z.size(); i++){
                std::cout << this->m_bit_vector_z[i] << ", ";
            }
            std::cout << std::endl;
        }

        iterator begin() const
        {
            return iterator(this, 0);
        }

        iterator end() const
        {
            return iterator(this, size());
        }

        size_type size()const{
            return m_size;
        }

        void swap(oz_vector &v) {

            m_bit_vector_o.swap(v.m_bit_vector_o);
            m_bit_vector_z.swap(v.m_bit_vector_z);
            sdsl::util::swap_support(m_support_select_o, v.m_support_select_o, &m_bit_vector_o, &v.m_bit_vector_o);
            sdsl::util::swap_support(m_support_select_z, v.m_support_select_z, &m_bit_vector_z, &v.m_bit_vector_z);
            sdsl::util::swap_support(m_support_rank_o, v.m_support_rank_o, &m_bit_vector_o, &v.m_bit_vector_o);
            sdsl::util::swap_support(m_support_rank_z, v.m_support_rank_z, &m_bit_vector_z, &v.m_bit_vector_z);
            std::swap(m_z_1, v.m_z_1);
            std::swap(m_size, v.m_size);
            std::swap(m_ones, v.m_ones);
            std::swap(m_zeroes, v.m_zeroes);
            std::swap(m_right_bsearch, v.m_right_bsearch);
        }

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="")const
        {
            sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += sdsl::write_member(m_size, out, child, "size");
            written_bytes += sdsl::write_member(m_z_1, out, child, "z_1");
            written_bytes += sdsl::write_member(m_ones, out, child, "ones");
            written_bytes += sdsl::write_member(m_zeroes, out, child, "zeroes");
            written_bytes += sdsl::write_member(m_right_bsearch, out, child, "right_bsearch");
            written_bytes += m_bit_vector_o.serialize(out, child, "bit_vector_o");
            written_bytes += m_bit_vector_z.serialize(out, child, "bit_vector_z");
            written_bytes += m_support_select_o.serialize(out, child, "select_o");
            written_bytes += m_support_select_z.serialize(out, child, "select_z");
            written_bytes += m_support_rank_o.serialize(out, child, "rank_o");
            written_bytes += m_support_rank_z.serialize(out, child, "rank_z");
            sdsl::structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream& in)
        {
            sdsl::read_member(m_size, in);
            sdsl::read_member(m_z_1, in);
            sdsl::read_member(m_ones, in);
            sdsl::read_member(m_zeroes, in);
            sdsl::read_member(m_right_bsearch, in);
            m_bit_vector_o.load(in);
            m_bit_vector_z.load(in);
            m_support_select_o.load(in, &m_bit_vector_o);
            m_support_select_z.load(in, &m_bit_vector_z);
            m_support_rank_o.load(in, &m_bit_vector_o);
            m_support_rank_z.load(in, &m_bit_vector_z);
        }

    };

    template<uint8_t t_b>
    class select_support_oz {

    public:
        typedef sdsl::bit_vector::size_type size_type;
        typedef oz_vector<> bit_vector_type;
    private:
        const bit_vector_type* m_v;
    public:

        explicit select_support_oz(const bit_vector_type* v=nullptr)
        {
            set_vector(v);
        }

        //! Returns the position of the i-th occurrence in the bit vector.
        size_type select(size_type i)const;


        size_type operator()(size_type i)const
        {
            return select(i);
        }

        size_type size()const
        {
            return m_v->size();
        }

        void set_vector(const bit_vector_type* v=nullptr)
        {
            m_v = v;
        }

        select_support_oz& operator=(const select_support_oz& ss)
        {
            if (this != &ss) {
                set_vector(ss.m_v);
            }
            return *this;
        }

        void swap(select_support_oz&) { }

        void load(std::istream&, const bit_vector_type* v=nullptr)
        {
            set_vector(v);
        }

        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="")const
        {
            return serialize_empty_object(out, v, name, this);
        }


    };

    template<uint8_t t_b>
    class succ_support_oz {

    public:
        typedef sdsl::bit_vector::size_type size_type;
        typedef oz_vector<> bit_vector_type;
    private:
        const bit_vector_type* m_v;

        size_type binary_search(size_type value)const{
            size_type l , r, m;
            l = 1; r = m_v->right_bsearch+1;
            size_type sum = 0;
            while (l<r) {
                m=(l+r)/2;
                sum = m_v->support_select_z(m) +m_v->support_select_o(m);
                if (sum < value)
                    l=m+1;
                else if (sum == value)
                    r=l=m+1;
                else
                    r=m;
            }
            return l-1;
        }

    public:

        explicit succ_support_oz(const bit_vector_type* v=nullptr)
        {
            set_vector(v);
        }

        size_type succ(size_type i)const;

        size_type operator()(size_type i)const
        {
            return succ(i);
        }

        size_type size()const
        {
            return m_v->size();
        }

        void set_vector(const bit_vector_type* v=nullptr)
        {
            m_v = v;
        }

        succ_support_oz& operator=(const succ_support_oz& ss)
        {
            if (this != &ss) {
                set_vector(ss.m_v);
            }
            return *this;
        }

        void swap(succ_support_oz&) {
        }

        void load(std::istream&, const bit_vector_type* v=nullptr)
        {
            set_vector(v);
        }

        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="")const
        {
            return serialize_empty_object(out, v, name, this);
        }


    };

    template<uint8_t t_b>
    class rank_support_oz {

    public:
        typedef sdsl::bit_vector::size_type size_type;
        typedef oz_vector<> bit_vector_type;
    private:
        const bit_vector_type* m_v;
    public:

        explicit rank_support_oz(const bit_vector_type* v=nullptr)
        {
            set_vector(v);
        }

        size_type binary_search(size_type value)const{
            size_type l , r, m;
            l = 1; r = m_v->right_bsearch+1;
            size_type sum = 0;
            while (l<r) {
                m=(l+r)/2;
                sum = m_v->support_select_z(m) +m_v->support_select_o(m);
                if (sum < value)
                    l=m+1;
                else if (sum == value)
                    r=l=m+1;
                else
                    r=m;
            }
            return l-1;
        }

        //! Returns the position of the i-th occurrence in the bit vector.
        size_type rank(size_type i)const;

        size_type operator()(size_type i)const
        {
            return rank(i-1);
        }

        size_type size()const
        {
            return m_v->size();
        }

        void set_vector(const bit_vector_type* v=nullptr)
        {
            m_v = v;
        }

        rank_support_oz& operator=(const rank_support_oz& ss)
        {
            if (this != &ss) {
                set_vector(ss.m_v);
            }
            return *this;
        }

        void swap(rank_support_oz&) { }

        void load(std::istream&, const bit_vector_type* v=nullptr)
        {
            set_vector(v);
        }

        size_type serialize(std::ostream& out, sdsl::structure_tree_node* v=nullptr, std::string name="")const
        {
            return serialize_empty_object(out, v, name, this);
        }


    };


    template<>
    typename rank_support_oz<0>::size_type rank_support_oz<0>::rank(size_type i)const
    {
        assert(i >= 0); assert(i <= m_v->size());
        if (i == 0) { return 0; }
        --i;
        size_type j = binary_search(i);
        if(m_v->z_1){
            if(i<(m_v->support_select_z(j+1) +m_v->support_select_o(j))){
                return i - m_v->support_select_o(j) + 1; //j-th run zeroes
            }else{
                return m_v->support_select_z(j+1); //j-th run ones
            }

        }else{
            if(i< (m_v->support_select_o(j+1) + m_v->support_select_z(j))){
                return m_v->support_select_z(j);
            }else{
                return i - m_v->support_select_o(j+1)+1;
            }
        }

    }



    template<>
    typename rank_support_oz<1>::size_type rank_support_oz<1>::rank(size_type i)const
    {
        assert(i >= 0); assert(i <= m_v->size());
        if (i == 0) { return 0; }
        --i;
        size_type j = binary_search(i);
        if(m_v->z_1){
            if(i<(m_v->support_select_z(j+1) + m_v->support_select_o(j))){
                return m_v->support_select_o(j); //j-th run zeroes
            }else{
                return i - m_v->support_select_z(j+1) + 1;//j-th run ones
            }

        }else{
            if(i< (m_v->support_select_o(j+1) + m_v->support_select_z(j))){
                return i - m_v->support_select_z(j)+1; //j-th run ones
            }else{
                return m_v->support_select_o(j+1);//j-th run zeroes
            }
        }
    }

    template<>
    typename select_support_oz<0>::size_type  select_support_oz<0>::select(size_type i)const
    {
        assert(i > 0);
        if(i > m_v->zeroes) return m_v->size();
        size_type rnk = m_v->support_rank_z(i);
        size_type sel = m_v->support_select_o(rnk+1-m_v->z_1);
        return i + sel -1;
    }

    template<>
    typename select_support_oz<1>::size_type  select_support_oz<1>::select(size_type i)const
    {
        assert(i > 0);
        if(i > m_v->ones) return m_v->size();
        size_type rnk = m_v->support_rank_o(i);
        size_type sel = m_v->support_select_z(rnk+m_v->z_1);
        return i + sel -1;
    }

    template<>
    typename succ_support_oz<1>::size_type succ_support_oz<1>::succ(size_type i) const {
        size_type j = binary_search(i);
        size_type k;
        if(m_v->z_1){
            k = m_v->support_select_z(j+1) +m_v->support_select_o(j);
            if(i< k){
                return k; //j-th run zeroes
            }else{
                return i; //j-th run ones
            }

        }else{
            k = m_v->support_select_o(j+1) + m_v->support_select_z(j);
            if(i< k){
                return i;  //j-th run ones
            }else{
                return k; //j-th run zeroes
            }
        }
    }


    template<>
    typename succ_support_oz<0>::size_type succ_support_oz<0>::succ(size_type i) const {
        size_type j = binary_search(i);
        size_type k;
        if(m_v->z_1){
            k = m_v->support_select_z(j+1) +m_v->support_select_o(j);
            if(i< k){
                return i; //j-th run zeroes
            }else{
                return k; //j-th run ones
            }

        }else{
            k = m_v->support_select_o(j+1) + m_v->support_select_z(j);
            if(i< k){
                return k;  //j-th run ones
            }else{
                return i; //j-th run zeroes
            }
        }
    }
}


#endif //RCT_OZ_VECTOR_HPP
