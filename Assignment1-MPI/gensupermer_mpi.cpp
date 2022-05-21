// *********************************************************************
//     [NAME]: QIU Yaowen      [STUDENT ID]: 20784389
//     [EMAIL]: yqiuau@connect.ust.hk
//     NOTICE: Write your code only in the specified section.
// *********************************************************************
// 7 MAR (update2.1), 28 FEB (update1): UPDATES IN read2supermers(...)
#define _in_
#define _out_
#define _MPI_TEST_
// #define DEBUG
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <chrono>
#include "utilities.hpp"
#ifdef _MPI_TEST_
#include "mpi.h"
#endif
using namespace std;

// void read2supermers(const char* _read, int read_len, int k, int p, _out_ char* &supermers, _out_ /**/int* &supermer_offs, int &n_supermers);
void read2supermers(const char *_read, int read_len, int k, int p, _out_ vector<string> &supermers);

const int MAX_PROCESS = 64;
int K, P;

int main(int argc, char **argvs)
{
#ifdef _MPI_TEST_
    MPI_Init(&argc, &argvs);
    MPI_Comm comm;
    int num_process; // number of processors
    int my_rank;     // my global rank
    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &num_process);
    MPI_Comm_rank(comm, &my_rank);
#endif

    int correctness_checking = 0;
    string output_path;
    string read_file_path;
    ArgParser(argc, argvs, K, P, read_file_path, correctness_checking, output_path);
    vector<string> reads;

    // Input data (the reads in CSR format)
    int num_of_reads = 0;
    char *reads_CSR;
    /**/ int *reads_CSR_offs;

    // Output data, each supermers should be a string in the vector
    // you need to save all the supermers to the vector below in the root(0) process
    vector<string> all_supermers;

#ifdef _MPI_TEST_
    if (my_rank == 0)
    {
        cout << "MPI running with " << num_process << " threads." << endl
             << endl;
#endif
        LoadReadsFromFile(read_file_path.c_str(), reads);
        Vector2CSR(reads, num_of_reads, reads_CSR, reads_CSR_offs);
        cout << reads.size() << " reads loaded from " << read_file_path << endl
             << endl;
#ifdef _MPI_TEST_
    }
#endif

    // time measurement starts
    auto start_time = chrono::high_resolution_clock::now();

#ifdef _MPI_TEST_
    //       hint: Input: "num_of_reads", "reads_CSR", "reads_CSR_offs", "K", "P"
    //       You need to save all the generated supermers in the vector "all_supermers" in Process 0.
    //       you need to do:
    //       1. Scatter the read data to each MPI processes.
    //       2. Perform the super-mer generation in each process.
    //          (you can refer to the sequential version to know the usage of the function read2supermers(...))
    //       3. Gather all the super-mers to the root process and store in the vector "all_supermers". (The order in the vector doesn't matter.)

    // ==============================================================
    // ==============================================================
    // ====    Write your implementation only below this line    ====
    // ==============================================================

    // Define Necessary variables for programming
    bool is_master = (my_rank == 0); // define master node signal
    int index;                       // initial value
    int bound;                       // upper bound value

    /*
    START SCATTER READS_EACH_PROCESSORS
    */
    int local_num_of_reads;     // local variable for each processor to store given counts
    int *reads_each_processors; // array to store nums of reads each processors
    if (is_master)
    {

        int mod = num_of_reads % num_process;      // one processor have mod num of reads only
        int quotient = num_of_reads / num_process; // num of reads per processor on average

        reads_each_processors = (int *)malloc(num_process * sizeof(int));

        index = 0;
        bound = num_process;
        while (index < bound)
        {

            if (index == 0)
            {
                reads_each_processors[index] = quotient + mod;
            }

            else
            {
                reads_each_processors[index] = quotient;
            }

            index++;
        }
    }

    // Scatter the array of reads_each_processors to local processes
    MPI_Scatter(reads_each_processors, 1, MPI_INT, &local_num_of_reads, 1, MPI_INT, 0, comm);

    /*
    END SCATTER READS_EACH_PROCESSORS
    */

    /*
    START SCATTER READS_CSR_OFFS
    */
    int *local_read_CSR_offs; // local variable to store read_CSR_Offs array
    int *nums_read_CSR_offs;  // array to indicate how many read_CSR_Offs for each processor
    int *displs;              // the array of pointers point to start of each read

    if (is_master)
    {
        displs = (int *)malloc(num_process * sizeof(int));             // allocate  array
        nums_read_CSR_offs = (int *)malloc(num_process * sizeof(int)); // allocate  array

        index = 0;
        bound = num_process;
        while (index < bound)
        {

            nums_read_CSR_offs[index] = reads_each_processors[index] + 1;

            if (index > 0)
            {
                displs[index] = displs[index - 1] + reads_each_processors[index - 1]; // the disls update formular
            }

            else
            {
                displs[index] = 0; // first element start from zero
            }
            index++;
        }
    }

    // since offset array have one extra pointer in the end, +1 is a must
    int size_local_read_CSR_offs = local_num_of_reads + 1;
    local_read_CSR_offs = (int *)malloc((size_local_read_CSR_offs) * sizeof(int));

    MPI_Scatterv(reads_CSR_offs, nums_read_CSR_offs, displs, MPI_INT, local_read_CSR_offs, size_local_read_CSR_offs, MPI_INT, 0, comm);

    /*
    END SCATTER READS_CSR_OFFS
    */

    /*
    START SCATTER READS_CSR
    */

    index = 0;
    bound = local_num_of_reads;
    bool flag = (reads.size() == 5);             // for testing
    int subtract_value = local_read_CSR_offs[0]; // record first value for futher processing

    while (index <= bound)
    {
        if (index == 0)
        {
            local_read_CSR_offs[index] = 0; //first must be zero
        }
        else
        {
            local_read_CSR_offs[index] = local_read_CSR_offs[index] - subtract_value; //perform left shift since spiting 
        }
        index++;
    }

    char *local_reads_CSR; // store scattered reads_CSR string vector
    int sizeof_local_reads_CSR = local_read_CSR_offs[local_num_of_reads] + 1; //+1 to keep poisiton for last pointer
    local_reads_CSR = (char *)malloc(sizeof_local_reads_CSR);

    int *local_sizeof_read_CSR;     // size of scattered reads_CSR string vector
    int *local_sizeof_readCSR_offs; // size of scattered readCSR_offs

    if (is_master)
    {
        local_sizeof_readCSR_offs = (int *)malloc(num_process * sizeof(int));

        local_sizeof_read_CSR = (int *)malloc(num_process * sizeof(int));

        // if (!flag)
        // {
        //     index = 1;
        // }
        // else
        // {
        index = 0;
        // }

        bound = num_process;
        while (index < bound)
        {
            if (index == 0)
            {
                local_sizeof_readCSR_offs[index] = 0;
            }
            else
            {
                local_sizeof_readCSR_offs[index] = reads_CSR_offs[displs[index]];

                // the size of read_CSR can be obtained from the gap value of a 2-sliding window
                local_sizeof_read_CSR[index - 1] = local_sizeof_readCSR_offs[index] - local_sizeof_readCSR_offs[index - 1];
            }

            index++;
        }

        local_sizeof_read_CSR[bound - 1] = reads_CSR_offs[num_of_reads] - local_sizeof_readCSR_offs[bound - 1];
    }

    MPI_Scatterv(reads_CSR, local_sizeof_read_CSR, local_sizeof_readCSR_offs, MPI_CHAR, local_reads_CSR, local_read_CSR_offs[local_num_of_reads], MPI_CHAR, 0, comm);

    /*
    END SCATTER READS_CSR
    */

    /*
    START generate supermers in local
    */

    vector<string> local_all_supermers;
    for (int index = 0; index < local_num_of_reads; index++)
    {
        vector<string> supermers_local;
        read2supermers(
            local_reads_CSR + local_read_CSR_offs[index],
            local_read_CSR_offs[index + 1] - local_read_CSR_offs[index],
            K, P,
            supermers_local);
        local_all_supermers.insert(local_all_supermers.end(), supermers_local.begin(), supermers_local.end());
    }

    /*
    END generate supermers in local
    */

    /*
    START GATHER number of local supermers in array to main process
    */

    char *local_supermers_CSR;
    int *local_supermers_CSR_OFFs;
    int nums_supermers;
    Vector2CSR(local_all_supermers, nums_supermers, local_supermers_CSR, local_supermers_CSR_OFFs);

    int *recvbuf;
    if (is_master)
    {
        recvbuf = (int *)malloc(sizeof(int) * num_process);
    }

    // Gather the count of supermers in each local process for futher use
    MPI_Gather(&nums_supermers, 1, MPI_INT, recvbuf, 1, MPI_INT, 0, comm);

    /*
    END GATHER number of local supermers in array to main process
    */

    /*
    START GATHER local supermers csroff to main process
    */
    int *local_supermers_CSR_OFFs_displs;  // total supermers offset displs
    int *local_supermers_CSR_OFFs_recvbuf; // total supermers offset receive buffer
    int *sizeof_local_supermers_CSR_OFFs;  // number of each local supermers in array

    int sizeof_total_local_supermers_CSR_OFFs = 0; // number of local supermers
    int sizeof_total_supermers = 0;                // number of aggreagred supermers in main process (should be)

    if (is_master)
    {
        local_supermers_CSR_OFFs_displs = (int *)malloc(sizeof(num_process) * num_process);
        sizeof_local_supermers_CSR_OFFs = (int *)malloc(sizeof(num_process) * num_process);

        index = 0;
        bound = num_process;
        while (index < bound)
        {
            sizeof_total_supermers = sizeof_total_supermers + recvbuf[index]; // just add to count total supers
            sizeof_local_supermers_CSR_OFFs[index] = recvbuf[index] + 1;      // offset points equals to nums of counter + 1

            if (index == 0)
            {
                local_supermers_CSR_OFFs_displs[index] = 0;
            }

            else
            {
                // by the definition of displs pointer
                local_supermers_CSR_OFFs_displs[index] = sizeof_local_supermers_CSR_OFFs[index - 1] + local_supermers_CSR_OFFs_displs[index - 1];
            }

            index++;
        }

        sizeof_total_local_supermers_CSR_OFFs = sizeof_total_supermers + num_process;

        local_supermers_CSR_OFFs_recvbuf = (int *)malloc((sizeof(int) * sizeof_total_local_supermers_CSR_OFFs));
    }

    MPI_Gatherv(local_supermers_CSR_OFFs, nums_supermers + 1, MPI_INT, local_supermers_CSR_OFFs_recvbuf, sizeof_local_supermers_CSR_OFFs, local_supermers_CSR_OFFs_displs, MPI_INT, 0, comm);

    /*
    END GATHER local supermers csroff to main process
    */

    /*
    START GATHER local supermers CSR to main process
    */

    int *total_supermers_nums;          // counter of all supermers
    int *total_supermers_lens;          // length of char of all supermers
    int *total_supermers_read_CSR_OFFs; // array to csr offset of all supermers

    if (is_master)
    {
        int recorder = 0; // tmp variable

        total_supermers_lens = (int *)malloc(sizeof_total_supermers * sizeof(int));

        total_supermers_read_CSR_OFFs = (int *)malloc((sizeof_total_supermers + 1) * sizeof(int));
        total_supermers_read_CSR_OFFs[0] = 0;

        total_supermers_nums = (int *)malloc(num_process * sizeof(int));

        displs = (int *)malloc(num_process * sizeof(int));
        displs[0] = 0; // first must be zero

        int offset_ptr = 0;          // index pointer for readcsr offs array
        int displs_ptr = offset_ptr; // index pointer for displs array

        index = 1;
        bound = sizeof_total_local_supermers_CSR_OFFs;
        while (index < bound)
        {
            if (local_supermers_CSR_OFFs_recvbuf[index] == 0) //
            {
                recorder = total_supermers_read_CSR_OFFs[offset_ptr]; // record value

                displs_ptr += 1;               // move to next index
                displs[displs_ptr] = recorder; // move to the index of next beginning sequence
                total_supermers_nums[displs_ptr - 1] = displs[displs_ptr] - displs[displs_ptr - 1];
            }
            else
            {
                offset_ptr += 1; // move to next index
                total_supermers_read_CSR_OFFs[offset_ptr] = local_supermers_CSR_OFFs_recvbuf[index] + recorder;
                total_supermers_lens[offset_ptr - 1] = total_supermers_read_CSR_OFFs[offset_ptr] - total_supermers_read_CSR_OFFs[offset_ptr - 1];
            }

            index++;
        }
        total_supermers_nums[displs_ptr] = total_supermers_read_CSR_OFFs[sizeof_total_supermers] - displs[num_process - 1];
    }

    char *recvbuffer; // reciever buffer to store aggregated supermers, follows the doc

    if (is_master)
    {
        // assign array for buffer
        recvbuffer = (char *)malloc(total_supermers_read_CSR_OFFs[sizeof_total_supermers] * sizeof(char));
    }

    MPI_Gatherv(local_supermers_CSR, local_supermers_CSR_OFFs[nums_supermers], MPI_CHAR, recvbuffer, total_supermers_nums, displs, MPI_CHAR, 0, comm);

    /*
    END GATHER local supermers CSR to main process
    */

    /*
    recvbuffer to all_supermers (vector string)
    */

    if (is_master)
    {

        string all_supermers_string = recvbuffer; // covert recvbuffer to string

        index = 0;
        bound = sizeof_total_supermers;

        while (index < bound)
        {

            // for each substr split by offset and len, append it into all_supermers vector
            // substr method can extract the sequence of string [begin index, end index]
            all_supermers.push_back(all_supermers_string.substr(total_supermers_read_CSR_OFFs[index], total_supermers_lens[index]));

            index++;
        }
    }

// ==============================================================
// ====    Write your implementation only above this line    ====
// ==============================================================
// ==============================================================
#endif

#ifdef _MPI_TEST_
    if (my_rank == 0)
    {
#endif
        // time measurement ends
        auto end_time = chrono::high_resolution_clock::now();
        auto duration_sec = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count() / 1000.0;
        cout << "Your algorithm finished in " << duration_sec << " sec." << endl
             << endl;

        // output to text file and correctness checking
        delete reads_CSR;
        delete reads_CSR_offs;
        if (correctness_checking)
            CorrectnessChecking(reads, K, P, all_supermers);
        if (!output_path.empty())
        {
            if (!correctness_checking)
                sort(all_supermers.begin(), all_supermers.end());
            SaveSupermers(output_path, all_supermers);
        }
#ifdef _MPI_TEST_
    }
    MPI_Barrier(comm);
    // cout<<"Thread "<<my_rank<<" ends."<<endl;

    MPI_Finalize();
#endif

    return 0;
}

/*
This function receives a C-style read string, the length of the read,
k (length of k-mer), p (length of miniizer), and output the supermers
which can be generated from this read to a vector<string>.
*/
void read2supermers(const char *_read, int read_len, int k, int p, _out_ vector<string> &supermers)
{
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
    for (i = p; i < k; i++)
    {
        new_minimizer = new_minimizer.substr(1, p - 1) + read[i]; // UPDATE1
        if (new_minimizer <= minimizer)
            minimizer = new_minimizer, mm_begin_pos = i - p + 1;
    }

    // Continue generating minimizers:
    for (i = 1; i < read_len - k + 1; i++)
    { // i: the beginning position of the current k-mer
        if (i > mm_begin_pos)
        {
            // new minimizer required
            prev_minimizer = minimizer;
            minimizer = new_minimizer = read.substr(i, p);
            for (j = i + p; j < i + k; j++)
            {
                new_minimizer = new_minimizer.substr(1, p - 1) + read[j]; // UPDATE1
                if (new_minimizer <= minimizer)
                    minimizer = new_minimizer, mm_begin_pos = j - p + 1;
            }
            // if the new minimizer equals to the previous one, we can continue
            if (minimizer != prev_minimizer)
            {
                skm_end_pos = i - 1 + k;
                supermers.push_back(read.substr(skm_begin_pos, skm_end_pos - skm_begin_pos)); // save the supermer
                skm_begin_pos = i;
            }
        }
        else
        {
            new_minimizer = read.substr(i + k - p, p); // UPDATE1
            if (new_minimizer < minimizer)
            { // save the supermer
                skm_end_pos = i - 1 + k;
                supermers.push_back(read.substr(skm_begin_pos, skm_end_pos - skm_begin_pos));
                skm_begin_pos = i;
                minimizer = new_minimizer, mm_begin_pos = i + k - 1 - p + 1;
            }
            if (new_minimizer == minimizer)
                mm_begin_pos = i + k - 1 - p + 1; // UPDATE1
        }
    } // UPDATE 2.1
    skm_end_pos = read_len;
    supermers.push_back(read.substr(skm_begin_pos, skm_end_pos - skm_begin_pos));
}