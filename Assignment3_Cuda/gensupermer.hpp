#pragma once

#include <cmath>
#include <cstring>
#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <string>
#include <cuda_runtime_api.h>
#include <cuda.h>

using namespace std;

inline void GPUAssert(cudaError_t code, const char *file, int line, bool abort = true) {
    if (code != cudaSuccess) {
        fprintf(stderr, "GPU assert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort)
            exit(code);
    }
};

#define GPUErrChk(ans) { GPUAssert((ans), __FILE__, __LINE__); }


void GenerateSupermer_GPU(vector<string> &reads, int K, int P, vector<string> &all_supermers, int num_blocks_per_grid, int num_threads_per_block);
