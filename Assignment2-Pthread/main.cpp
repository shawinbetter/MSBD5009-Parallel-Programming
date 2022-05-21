/* Keep this file unchanged
 * COMPILE: g++ -std=c++11 -lpthread gensupermer_pthread.cpp main.cpp -o pthread_gs
 *                  ps, in some platform, the command "-lpthread" should be replaced by "-pthread"
 * RUN:     ./pthread_gs <K> <P> <num_threads> <data_file> <correctness_check> <result_folder (optional)>
 */

#include <cassert>
#include <chrono>
#include <iostream>

#include "gensupermer.hpp"
#include "utilities.hpp"
using namespace std;


int main(int argc, char **argvs) {
    int K, P;
    int num_of_threads = 0;
    int correctness_checking = 0;
    string output_path;
    string read_file_path;
    if (argc < 6) {
        cerr << "Error! usage:" << endl << "\"./gs [k] [p] [num_of_threads] [data file path] [correctness_checking=0/1] (optional) [supermers_output_path]\"" << endl;
        exit(1);
    }
    K = atoi(argvs[1]);
    P = atoi(argvs[2]);
    num_of_threads = atoi(argvs[3]);
    read_file_path = argvs[4];
    if (argc >= 6) correctness_checking = atoi(argvs[5]);
    if (argc >= 7) {
        output_path = string(argvs[6]);
        if (*(output_path.end()-1) != '/') output_path += "/my_gs_output.txt";
        else output_path += "my_gs_output.txt";
    }
    vector<string> reads;
    
    // Input data (the reads in CSR format)
    int num_of_reads = 0;
    char* reads_CSR;
    /**/int* reads_CSR_offs;

    // Output data, each supermers should be a string in the vector
    vector<string> all_supermers;

    
    LoadReadsFromFile(read_file_path.c_str(), reads);
    Vector2CSR(reads, num_of_reads, reads_CSR, reads_CSR_offs);
    cout << reads.size() << " reads loaded from "<< read_file_path << endl << endl;
    
    // time measurement starts
    auto start_time = chrono::high_resolution_clock::now();
    gensupermers(reads_CSR, reads_CSR_offs, K, P, num_of_reads, all_supermers, num_of_threads);
    // time measurement ends
    auto end_time = chrono::high_resolution_clock::now();
    auto duration_sec = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count()/1000.0;
    printf("Your algorithm finished in %.5g sec.\n", duration_sec);
    
    // output to text file and correctness checking
    delete reads_CSR;
    delete reads_CSR_offs;
    if (correctness_checking) CorrectnessChecking(reads, K, P, all_supermers);
    if (!output_path.empty()) {
        if (!correctness_checking) sort(all_supermers.begin(), all_supermers.end());
        SaveSupermers(output_path, all_supermers);
    }
    
    return 0;
}
