#!/bin/bash

# Ensure all commands succeed.
set -e

# Get start timestamp.
start_datetime=$(date +'%Y%m%d_%H%M%S')

# Define directories.
output_dir=out
alternatives_dir=alternatives

# Create directories.
mkdir -p ${output_dir}
rm -f ${output_dir}/*

# Define files.
experiment_out_plcp=${output_dir}/experiment_plcp.out
experiment_out_bwt_select=${output_dir}/experiment_bwt_select.out

# Clear output and version file.
rm -f ${experiment_out_plcp}
rm -f ${experiment_out_bwt_select}

# Initialize submodules and files.
make setup

# Apply fixes.
cp alternatives/external/zombit-vector/include/oz_vector__00000000_000000.hpp external/zombit-vector/include/oz_vector.hpp
cp alternatives/external/zombit-vector/include/partitioned_zombit_vector_sparse__00000000_000000.hpp external/zombit-vector/include/partitioned_zombit_vector_sparse.hpp

# Run test
index_type_list=("sdsl::hyb_vector<8>" "sdsl::hyb_vector<16>" "sdsl::hyb_vector<32>" "sdsl::hyb_vector<64>")
index_alias_list=("HYB" "HYB" "HYB" "HYB")
block_size_list=("8" "16" "32" "64")
additional_flags_list=("" "" "" "")
for ((index_type_index_alias_block_size_additional_flags_i=0; index_type_index_alias_block_size_additional_flags_i<4; ++index_type_index_alias_block_size_additional_flags_i)); do
index_type=${index_type_list[$index_type_index_alias_block_size_additional_flags_i]}
index_alias=${index_alias_list[$index_type_index_alias_block_size_additional_flags_i]}
block_size=${block_size_list[$index_type_index_alias_block_size_additional_flags_i]}
additional_flags=${additional_flags_list[$index_type_index_alias_block_size_additional_flags_i]}
make  CPPFLAGS="-DTEST_LARGE_FILES -DINDEX_TYPE=\"${index_type}\" -DINDEX_ALIAS=\"${index_alias}\" -DBLOCK_SIZE=\"${block_size}\" ${additional_flags}" clean build/test_select
./build/test_select
if [ ${PIPESTATUS[0]} -ne 0 ]; then exit 1; fi
done

# Print header
make  CPPFLAGS="-DNDEBUG -DHEADERLESS -DFOOTERLESS" clean build/bench_header
./build/bench_header | tee -a "${experiment_out_plcp}"
if [ ${PIPESTATUS[0]} -ne 0 ]; then exit 1; fi

# PLCP Benchmark
for input_file in "files/dna.200MB" "files/proteins.200MB" "files/english.200MB" "files/sources.200MB" "files/dblp.xml.200MB" "files/para" "files/cere" "files/influenza" "files/world_leaders" "files/kernel"; do
index_type_list=("sdsl::rrr_vector<15>" "sdsl::rrr_vector<31>" "sdsl::rrr_vector<63>" "sdsl::rrr_vector<127>" "sdsl::hyb_vector<8>" "sdsl::hyb_vector<16>" "sdsl::hyb_vector<32>" "sdsl::hyb_vector<64>" "sdsl::sd_vector<>" "sdsl::bit_vector" "runs_vectors::oz_vector<>" "ZombitVectorWrapper<>" "NaiveHybVectorWrapper<sdsl::hyb_vector<8>>" "NaiveHybVectorWrapper<sdsl::hyb_vector<16>>" "NaiveHybVectorWrapper<sdsl::hyb_vector<32>>" "NaiveHybVectorWrapper<sdsl::hyb_vector<64>>" "PastaVectorWrapper<>" "LAVectorWrapper<>")
index_alias_list=("RRR" "RRR" "RRR" "RRR" "HYB" "HYB" "HYB" "HYB" "SD" "BV" "OZ" "ZOM" "NH" "NH" "NH" "NH" "PASTA" "LA")
block_size_list=("15" "31" "63" "127" "8" "16" "32" "64" "" "" "" "" "8" "16" "32" "64" "" "")
additional_flags_list=("" "" "" "" "" "" "" "" "" "" "-DUSE_OZ -I\"external/zombit-vector/include\" -I\"build/zombit-vector/include\"" "-DUSE_ZOMBIT -I\"external/zombit-vector/include\" -I\"build/zombit-vector/include\"" "-DUSE_NAIVE_HYB -I\"external/zombit-vector/include\" -I\"build/zombit-vector/include\"" "-DUSE_NAIVE_HYB -I\"external/zombit-vector/include\" -I\"build/zombit-vector/include\"" "-DUSE_NAIVE_HYB -I\"external/zombit-vector/include\" -I\"build/zombit-vector/include\"" "-DUSE_NAIVE_HYB -I\"external/zombit-vector/include\" -I\"build/zombit-vector/include\"" "-DUSE_PASTA -I\"external/pasta/include\" -I\"external/pasta_utils/include\" -I\"external/tlx\"" "-DUSE_LA_VECTOR -I\"external/la_vector/include\" -pthread -fopenmp")
for ((index_type_index_alias_block_size_additional_flags_i=0; index_type_index_alias_block_size_additional_flags_i<18; ++index_type_index_alias_block_size_additional_flags_i)); do
index_type=${index_type_list[$index_type_index_alias_block_size_additional_flags_i]}
index_alias=${index_alias_list[$index_type_index_alias_block_size_additional_flags_i]}
block_size=${block_size_list[$index_type_index_alias_block_size_additional_flags_i]}
additional_flags=${additional_flags_list[$index_type_index_alias_block_size_additional_flags_i]}
make  CPPFLAGS="-DNDEBUG -DHEADERLESS -DFOOTERLESS -DSEED=0 -DINDEX_TYPE=\"${index_type}\" -DINDEX_ALIAS=\"${index_alias}\" -DBLOCK_SIZE=\"${block_size}\" ${additional_flags}" clean build/bench_plcp
./build/bench_plcp ${input_file} | tee -a "${experiment_out_plcp}"
if [ ${PIPESTATUS[0]} -ne 0 ]; then exit 1; fi
done
done

# Print footer
make  CPPFLAGS="-DNDEBUG -DHEADERLESS -DFOOTERLESS" clean build/bench_footer
./build/bench_footer | tee -a "${experiment_out_plcp}"
if [ ${PIPESTATUS[0]} -ne 0 ]; then exit 1; fi

# Print header
make  CPPFLAGS="-DNDEBUG -DHEADERLESS -DFOOTERLESS" clean build/bench_header
./build/bench_header | tee -a "${experiment_out_bwt_select}"
if [ ${PIPESTATUS[0]} -ne 0 ]; then exit 1; fi

# BWT Select Benchmark
for input_file in "files/dna.200MB" "files/proteins.200MB" "files/english.200MB" "files/sources.200MB" "files/dblp.xml.200MB" "files/para" "files/cere" "files/influenza" "files/world_leaders" "files/kernel"; do
index_type_list=("sdsl::rrr_vector<15>" "sdsl::rrr_vector<31>" "sdsl::rrr_vector<63>" "sdsl::rrr_vector<127>" "sdsl::hyb_vector<8>" "sdsl::hyb_vector<16>" "sdsl::hyb_vector<32>" "sdsl::hyb_vector<64>" "sdsl::sd_vector<>" "sdsl::bit_vector" "runs_vectors::oz_vector<>" "ZombitVectorWrapper<>" "NaiveHybVectorWrapper<sdsl::hyb_vector<8>>" "NaiveHybVectorWrapper<sdsl::hyb_vector<16>>" "NaiveHybVectorWrapper<sdsl::hyb_vector<32>>" "NaiveHybVectorWrapper<sdsl::hyb_vector<64>>" "PastaVectorWrapper<>")
index_alias_list=("RRR" "RRR" "RRR" "RRR" "HYB" "HYB" "HYB" "HYB" "SD" "BV" "OZ" "ZOM" "NH" "NH" "NH" "NH" "PASTA")
block_size_list=("15" "31" "63" "127" "8" "16" "32" "64" "" "" "" "" "8" "16" "32" "64" "")
additional_flags_list=("" "" "" "" "" "" "" "" "" "" "-DUSE_OZ -I\"external/zombit-vector/include\" -I\"build/zombit-vector/include\"" "-DUSE_ZOMBIT -I\"external/zombit-vector/include\" -I\"build/zombit-vector/include\"" "-DUSE_NAIVE_HYB -I\"external/zombit-vector/include\" -I\"build/zombit-vector/include\"" "-DUSE_NAIVE_HYB -I\"external/zombit-vector/include\" -I\"build/zombit-vector/include\"" "-DUSE_NAIVE_HYB -I\"external/zombit-vector/include\" -I\"build/zombit-vector/include\"" "-DUSE_NAIVE_HYB -I\"external/zombit-vector/include\" -I\"build/zombit-vector/include\"" "-DUSE_PASTA -I\"external/pasta/include\" -I\"external/pasta_utils/include\" -I\"external/tlx\"")
for ((index_type_index_alias_block_size_additional_flags_i=0; index_type_index_alias_block_size_additional_flags_i<17; ++index_type_index_alias_block_size_additional_flags_i)); do
index_type=${index_type_list[$index_type_index_alias_block_size_additional_flags_i]}
index_alias=${index_alias_list[$index_type_index_alias_block_size_additional_flags_i]}
block_size=${block_size_list[$index_type_index_alias_block_size_additional_flags_i]}
additional_flags=${additional_flags_list[$index_type_index_alias_block_size_additional_flags_i]}
wt_type_list=("sdsl::wt_huff" "sdsl::wt_blcd")
wt_alias_list=("HUFF" "BLCD")
for ((wt_type_wt_alias_i=0; wt_type_wt_alias_i<2; ++wt_type_wt_alias_i)); do
wt_type=${wt_type_list[$wt_type_wt_alias_i]}
wt_alias=${wt_alias_list[$wt_type_wt_alias_i]}
make  CPPFLAGS="-DNDEBUG -DHEADERLESS -DFOOTERLESS -DSEED=0 -DINDEX_TYPE=\"${wt_type}<${index_type}>\" -DINDEX_ALIAS=\"${wt_alias}_${index_alias}\" -DBLOCK_SIZE=\"${block_size}\" ${additional_flags}" clean build/bench_bwt_select
./build/bench_bwt_select ${input_file} | tee -a "${experiment_out_bwt_select}"
if [ ${PIPESTATUS[0]} -ne 0 ]; then exit 1; fi
done
done
done

# Print footer
make  CPPFLAGS="-DNDEBUG -DHEADERLESS -DFOOTERLESS" clean build/bench_footer
./build/bench_footer | tee -a "${experiment_out_bwt_select}"
if [ ${PIPESTATUS[0]} -ne 0 ]; then exit 1; fi

# Get end timestamp.
end_datetime=$(date +'%Y%m%d_%H%M%S')

# Print timestamps and output location.
echo "Start time: ${start_datetime}"
echo "End time: ${end_datetime}"
echo "Output saved to ${output_dir}"

