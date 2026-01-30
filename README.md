# Prerequisites
- make
- cmake (to install SDSL)
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
- Load all necessary submodules.
- Install necessary libraries in `build`.
- Build and run `build/test_select`.
- Build and run `build/bench_plcp`, with the output written to `out/experiment_plcp.out`.
- Build and run `build/bench_bwt_select`, with the output written to `out/experiment_bwt_select.out`.
