// *********************************************************************
//     NOTICE: This is the sequential-version for your reference.
//             You can know the usage of the function "read2supermers"
//             from this file.
//     # Compile: g++ -std=c++11 gensupermer_sequential.cpp -o seq_gs
//     # ./seq_gs <K> <P> <data_file> <result_folder (optional)>
// *********************************************************************

#define _in_
#define _out_
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <chrono>
#include "utilities.hpp"
using namespace std;

// *********************************************************************
//     NOTICE: This is the sequential-version for your reference.
//             You can know the usage of the function "read2supermers"
//             from this file.
// *********************************************************************
// 7 MAR (update2.1), 28 FEB (update1): UPDATES IN read2supermers(...)
#define _in_
#define _out_
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <chrono>
#include "utilities.hpp"
using namespace std;

void read2supermers(const char* _read, int read_len, int k, int p, _out_ vector<string> &supermers);

const int MAX_PROCESS = 64;
int K, P;

int main(int argc, char **argvs) {

    string output_path;
    if (argc < 4) {
        cerr << "Error! usage:" << endl << "\"./seq_gs [k] [p] [data file path] (optional)[supermers_output_folder]\"" << endl;
        exit(1);
    }
    K = atoi(argvs[1]);
    P = atoi(argvs[2]);
    if (argc >= 5) output_path = string(argvs[4]) + string("seq_gs_output.txt");

    string read_file_path = string(argvs[3]);
    
    vector<string> reads;
    
    // Input data
    int num_of_reads = 0;
    char* reads_CSR;
    /**/int* reads_CSR_offs;

    // Output data, save all the supermers to the vector below in the root(0) process
    vector<string> all_supermers;

    LoadReadsFromFile(read_file_path.c_str(), reads);
    Vector2CSR(reads, num_of_reads, reads_CSR, reads_CSR_offs);
    cout << reads.size() << " reads loaded from "<< read_file_path << endl << endl;
    
    // time measurement starts
    auto start_time = chrono::high_resolution_clock::now();

    // ========== Reference ==========
    // Below is a sample of using the function "read2supermers" for your reference.
    #ifndef _MPI_TEST_
    for (int i=0; i<num_of_reads; i++) {
        vector<string> supermers_local;
        int new_supermer_size=0;
        read2supermers(
            reads_CSR + reads_CSR_offs[i],
            reads_CSR_offs[i+1] - reads_CSR_offs[i],
            K, P,
            supermers_local
        );
        all_supermers.insert(all_supermers.end(), supermers_local.begin(), supermers_local.end());
    }
    #endif
    // # Reference End #
    
    // time measurement ends
    auto end_time = chrono::high_resolution_clock::now();
    auto duration_sec = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count()/1000.0;
    cout << all_supermers.size() << " supermers were generated." << endl;
    cout << "Sequential algorithm finished in " << duration_sec << " sec." << endl << endl;
    
    // output to text file and correctness checking
    delete reads_CSR;
    delete reads_CSR_offs;
    if (!output_path.empty()) {
        sort(all_supermers.begin(), all_supermers.end());
        SaveSupermers(output_path, all_supermers);
    }
    
    return 0;
}

/*
This function receives a C-style read, the length of the read, k (length of k-mer), p (length of miniizer),
and output the supermers which can be generated from this read.
*/
void read2supermers(const char* _read, int read_len, int k, int p, _out_ vector<string> &supermers) {
    string prev_minimizer, minimizer, new_minimizer;
    string read(_read, read_len); // from-buffer init
    int i, j;
    char base;
    int skm_begin_pos, skm_end_pos, mm_begin_pos;
    
    // Generate the first k-mer's minimizer:
    skm_begin_pos = 0;
    skm_end_pos = k;
    mm_begin_pos = 0;
    minimizer = new_minimizer = read.substr(0, p);
    for (i=p; i<k; i++) {
        new_minimizer = new_minimizer.substr(1, p-1) + read[i]; // UPDATE1
        if (new_minimizer <= minimizer) minimizer = new_minimizer, mm_begin_pos = i-p+1;
    }

    // Continue generating minimizers:
    for (i=1; i<read_len-k+1; i++) { // i: the beginning position of the current k-mer
        if (i > mm_begin_pos) {
            // new minimizer required
            prev_minimizer = minimizer;
            minimizer = new_minimizer = read.substr(i, p);
            for (j=i+p; j<i+k; j++) {
                new_minimizer = new_minimizer.substr(1, p-1) + read[j]; // UPDATE1
                if (new_minimizer <= minimizer) minimizer = new_minimizer, mm_begin_pos = j-p+1;
            }
            // if the new minimizer equals to the previous one, we can continue
            if (minimizer != prev_minimizer) {
                skm_end_pos = i-1+k;
                supermers.push_back(read.substr(skm_begin_pos, skm_end_pos-skm_begin_pos)); // save the supermer
                skm_begin_pos = i;
            }
        }
        else {
            new_minimizer = read.substr(i+k-p, p); // UPDATE1
            if (new_minimizer < minimizer) { // save the supermer
                skm_end_pos = i-1+k;
                supermers.push_back(read.substr(skm_begin_pos, skm_end_pos-skm_begin_pos));
                skm_begin_pos = i;
                minimizer = new_minimizer, mm_begin_pos = i+k-1-p+1;
            }
            if (new_minimizer == minimizer) mm_begin_pos = i+k-1-p+1; // UPDATE1
        }
    } // UPDATE 2.1
    skm_end_pos = read_len;
    supermers.push_back(read.substr(skm_begin_pos, skm_end_pos-skm_begin_pos));
}
