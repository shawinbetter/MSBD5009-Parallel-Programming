/* Keep this file unchanged
 * COMPILE: nvcc -std=c++11 gensupermer_cuda.cu main.cpp  -o cuda
 * RUN:     ./cuda  <K> <P> <num_blocks_per_grid> <num_threads_per_block> <data_file> <correctness_check>
 */

#include <cassert>
#include <chrono>
#include <iostream>
#include "gensupermer.hpp"
#include "utilities.hpp"

using namespace std;

int main(int argc, char **argv) {
    int K, P;
    int num_blocks_per_grid = 0;
    int num_threads_per_block = 0;
    int correctness_checking = 0;
    string output_path;
    string read_file_path;
    
    
    if (argc < 7) {
        cerr << "Error! usage:" << endl << "\"./gs [k] [p] [num_blocks_per_grid] [num_threads_per_block]  [data file path] [correctness_checking=0/1] (optional) [supermers_output_path]\"" << endl;
        exit(1);
    }
    K = atoi(argv[1]);
    P = atoi(argv[2]);
    num_blocks_per_grid = atoi(argv[3]);
    num_threads_per_block = atoi(argv[4]);
    read_file_path = argv[5];
    if (argc >= 7) correctness_checking = atoi(argv[6]);
    if (argc >= 8) {
        output_path = string(argv[7]);
        if (*(output_path.end()-1) != '/') output_path += "/my_gs_output.txt";
        else output_path += "my_gs_output.txt";
    }
    vector<string> reads;
    

    // Output data, each supermers should be a string in the vector
    vector<string> all_supermers;

    
    LoadReadsFromFile(read_file_path.c_str(), reads);
    cout << reads.size() << " reads lala loaded from "<< read_file_path << endl << endl;


    cudaDeviceReset();
    // cudaEvent_t cuda_start, cuda_end;

    float kernel_time;
    auto start_clock = std::chrono::high_resolution_clock::now();
    
    // cudaEventCreate(&cuda_start);
    // cudaEventCreate(&cuda_end);

    // cudaEventRecord(cuda_start);

    GenerateSupermer_GPU(
        reads, K, P, all_supermers, num_blocks_per_grid, num_threads_per_block // update 1
    );

    // cudaEventRecord(cuda_end);

    // cudaEventSynchronize(cuda_start);
    // cudaEventSynchronize(cuda_end);

    // cudaEventElapsedTime(&kernel_time, cuda_start, cuda_end);
    // GPUErrChk(cudaDeviceSynchronize());

    auto end_clock = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> diff = end_clock - start_clock;
    
    printf("Elapsed Time: %.9lf s\n", diff.count());
    // fprintf(stderr, "Driver Time: %.9lf s\n", kernel_time / pow(10, 3));



    if (correctness_checking) CorrectnessChecking(reads, K, P, all_supermers);
    if (!output_path.empty()) {
        if (!correctness_checking) sort(all_supermers.begin(), all_supermers.end());
        SaveSupermers(output_path, all_supermers);
    }

    return 0;
}
