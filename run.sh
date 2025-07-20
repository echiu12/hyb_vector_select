make download
make tests/test_select && ./tests/test_select > tests/test_select.out
make tests/bench_plcp && ./tests/bench_plcp > tests/bench_plcp.out
make tests/bench_bwt_select && ./tests/bench_bwt_select > ./tests/bench_bwt_select
Rscript tests/visualize.R tests