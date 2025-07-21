# Prerequisites
- make
- cmake (to install SDSL)
- Rscript (with packages ggplot2, RColorBrewer, grid, and cowplot)
- wget (to get text files)
- gunzip (to get text files)
- g++

# Instructions
To run tests and benchmarks, run:
```
./run.sh
```
The script will:
- Download text files into `files`.
- Install SDSL into `sdsl`.
- Build and run `tests/test_select`, with the output written to `tests/test_select.out`.
- Build and run `tests/bench_plcp`, with the output written to `tests/bench_plcp.out`.
- Build and run `tests/bench_bwt_select`, with the output written to `tests/bench_bwt_select.out`.
- Create graphs for each `test/*.out` file written to `tests/*.pdf`.
